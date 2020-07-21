/* format_mapped_reads: a program to ensure SAM and BAM format reads
 * are conforming to expectations of methpipe software
 *
 * Copyright (C) 2020 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "cigar_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::runtime_error;
using std::unordered_map;
using std::swap;

static bool
is_mapped(const sam_rec &aln) {
  return true;
}

static bool
is_primary(const sam_rec &aln) {
  return true;
}

static bool
is_mapped_single_end(const sam_rec &aln) {
  return true;
}

static bool
is_t_rich(const sam_rec &aln) {
  return true;
}

static bool
is_a_rich(const sam_rec &aln) {
  return !is_t_rich(aln);
}


static void
flip_strand(sam_rec &aln) {
  if (check_flag(aln, samflags::read_rc))
    unset_flag(aln, samflags::read_rc);
  else
    set_flag(aln, samflags::read_rc);
}


static void
revcomp(sam_rec &aln) {
  flip_strand(aln); // set strand to opposite of current value
  revcomp_inplace(aln.seq); // reverse complement sequence
  std::reverse(begin(aln.qual), end(aln.qual)); // and quality scores
}

static uint32_t
get_r_end(const sam_rec &sr) {
  return sr.pos + cigar_rseq_ops(sr.cigar);
}

/********Below are functions for merging pair-end reads********/
static void
fill_overlap(const bool rc, const sam_rec &sr, const size_t start,
             const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = rc ? (get_r_end(sr) - end) : (start - sr.pos);
  const size_t b = rc ? (get_r_end(sr) - start) : (end -  sr.pos);
  copy(begin(sr.seq) + a, begin(sr.seq) + b, begin(seq) + offset);
  copy(begin(sr.qual) + a, begin(sr.qual) + b, begin(scr) + offset);
}

static size_t
merge_mates(const size_t suffix_len, const size_t range,
            const sam_rec &one, const sam_rec &two, sam_rec &merged) {

  const bool rc = check_flag(one, samflags::read_rc);
  const uint32_t overlap_left = std::max(one.pos, two.pos);
  const uint32_t overlap_right = min(get_r_end(one), get_r_end(two));

  uint32_t one_left = 0, one_right = 0, two_left = 0, two_right = 0;
  if (rc) {
    one_left = max(overlap_right, one.pos);
    one_right = get_r_end(one);

    two_left = two.pos;
    two_right = min(overlap_left, get_r_end(two));
  }
  else {
    one_left = one.pos;
    one_right = min(overlap_left, get_r_end(one));

    two_left = max(overlap_right, two.pos);
    two_right = get_r_end(two);
  }

  const int frag_len = rc ? (one_right - two_left) : (two_right - one_left);

  // assert(len > 0);
  // if the above assertion fails, it usually means the mair is
  // discordant (end1 downstream of end2). Also it means the SAM flag
  // of this pair of reads is not properly set. To avoid termination,
  // currently this assertion is ignored but no output will be
  // generated for discordant pairs.

  if (frag_len > 0) {
    merged.seq = string(frag_len, 'N');
    merged.qual = string(frag_len, 'B');
    if (frag_len <= static_cast<int>(range)) {
      // lim_one: offset in merged sequence where overlap starts
      const size_t lim_one = one_right - one_left;
      copy(begin(one.seq), begin(one.seq) + lim_one, begin(merged.seq));
      copy(begin(one.qual), begin(one.qual) + lim_one, begin(merged.qual));

      const size_t lim_two = two_right - two_left;
      copy(end(two.seq) - lim_two, end(two.seq), end(merged.seq) - lim_two);
      copy(end(two.qual) - lim_two, end(two.qual), end(merged.qual) - lim_two);

      // deal with overlapping part
      if (overlap_left < overlap_right) {
        const size_t one_bads = count(begin(one.seq), end(one.seq), 'N');
        const int info_one = one.seq.length() - (one_bads + one.mapq);

        const size_t two_bads = count(begin(two.seq), end(two.seq), 'N');
        const int info_two = two.seq.length() - (two_bads + two.mapq);

        // use the mate with the most info to fill in the overlap
        if (info_one >= info_two)
          fill_overlap(rc, one, overlap_left, overlap_right, lim_one,
                       merged.seq, merged.qual);
        else
          fill_overlap(rc, two, overlap_left, overlap_right, lim_one,
                       merged.seq, merged.qual);
      }
    }

    merged.pos = rc ? two.pos : one.pos;
    // merged.r.set_end(merged.pos + frag_len);
    const string name(one.qname);
    merged.qname = "FRAG:" + name.substr(0, name.size() - suffix_len);
  }
  return static_cast<size_t>(max(frag_len, 0));
}

inline static bool
is_same_read(const size_t suffix_len, const sam_rec &a, const sam_rec &b) {
  return std::equal(begin(a.qname), end(a.qname) - suffix_len, begin(b.qname));
}

/********Above are functions for merging pair-end reads********/

static string
remove_suffix(const sam_rec &aln, const size_t suffix_len) {
  return aln.qname.substr(0, aln.qname.size() - suffix_len);
}

static bool
precedes_by_more_than(const sam_rec &a, const sam_rec &b,
                      const size_t max_frag_len) {
  return (a.rname < b.rname ||
          (a.rname == b.rname && (get_r_end(a) + max_frag_len < b.pos)));
}


// ////////////////////////////////////////
// // BSMAP
// ////////////////////////////////////////

// class BSMAPFLAG : public FLAG {
// public:
//   BSMAPFLAG(const size_t f) : FLAG(f) {}
// };

// inline static void
// bsmap_get_strand(const string &strand_str, string &strand, string &bs_forward) {
//   strand = strand_str.substr(5, 1);
//   bs_forward = strand_str.substr(6, 1);
//   if (bs_forward == "-")
//     strand = (strand == "+" ? "-" : "+");
// }

// bool
// SAMReader::get_SAMRecord_bsmap(const string &str, sam_rec &samr) {

//   cerr << "WARNING: "<< "[BSMAP Converter] test "
//        << "version: may contain bugs" << endl;

//   string name, chrom, CIGAR, mate_name, seq, qual, strand_str, mismatch_str;
//   size_t flag, start, mapq_score, mate_start;
//   int seg_len;

//   std::istringstream iss(str);
//   if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
//         >> mate_name >> mate_start >> seg_len >> seq >> qual
//         >> mismatch_str >> strand_str)) {
//     good = false;
//     throw runtime_error("malformed line in bsmap SAM format:\n" + str);
//   }

//   BSMAPFLAG Flag(flag);

//   samr.mr.r.set_chrom(chrom);
//   samr.mr.r.set_start(start - 1);
//   samr.mr.r.set_name(name);
//   samr.mr.r.set_score(atoi(mismatch_str.substr(5).c_str()));

//   string strand, bs_forward;
//   bsmap_get_strand(strand_str, strand, bs_forward);
//   samr.mr.r.set_strand(strand[0]);

//   string new_seq, new_qual;
//   apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

//   samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size());
//   samr.mr.seq = new_seq;
//   samr.mr.qual = new_qual;

//   samr.is_Trich = Flag.is_Trich();
//   samr.is_mapping_paired = Flag.is_mapping_paired();

//   return good;
// }

////////////////////////////////////////
// Bismark
// ////////////////////////////////////////

// class BISMARKFLAG : public FLAG {
// public:
//   BISMARKFLAG(const size_t f) : FLAG(f) {}
//   bool is_Trich() const {return is_pairend() ? FLAG::is_Trich() : true;}
// };

// static size_t
// get_mismatch_bismark(const string &edit_distance_str,
//                      const string &meth_call_str) {
//   /* the result of this function might not be accurate, because if a
//   sequencing error occurs on a cytosine, then it probably will be
//   reported as a convertion
//   */
//   size_t edit_distance;
//   edit_distance = atoi(edit_distance_str.substr(5).c_str());

//   int convert_count = 0;
//   const char *temp = meth_call_str.substr(5).c_str();
//   while (*temp != '\0') {
//     if (*temp == 'x' || *temp == 'h' || *temp == 'z')
//       ++convert_count;
//     ++temp;
//   }

//   return edit_distance - convert_count;
// }

// bool
// SAMReader::get_SAMRecord_bismark(const string &str, sam_rec &samr) {
//   string name, chrom, CIGAR, mate_name, seq, qual, strand_str,
//     edit_distance_str, mismatch_str, meth_call_str,
//     read_conv_str, genome_conv_str;
//   size_t flag, start, mapq_score, mate_start;
//   int seg_len;

//   std::istringstream iss(str);
//   if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
//         >> mate_name >> mate_start >> seg_len >> seq >> qual
//         >> edit_distance_str >> mismatch_str >> meth_call_str
//         >> read_conv_str >> genome_conv_str)) {
//     good = false;
//     throw runtime_error("malformed line in bismark SAM format:\n" + str);
//   }

//   BISMARKFLAG Flag(flag);

//   samr.mr.r.set_chrom(chrom);
//   samr.mr.r.set_start(start - 1);
//   samr.mr.r.set_name(name);
//   samr.mr.r.set_score(get_mismatch_bismark(edit_distance_str, meth_call_str));
//   samr.mr.r.set_strand(Flag.is_revcomp() ? '-' : '+');

//   string new_seq, new_qual;
//   apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

//   if (Flag.is_revcomp()) {
//     revcomp_inplace(new_seq);
//     std::reverse(new_qual.begin(), new_qual.end());
//   }

//   samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size());
//   samr.mr.seq = new_seq;
//   samr.mr.qual = new_qual;

//   samr.is_Trich = Flag.is_Trich();
//   samr.is_mapping_paired = Flag.is_mapping_paired();

//   return good;
// }

// ////////////////////////////////////////
// // BS Seeker
// ////////////////////////////////////////

// class BSSEEKERFLAG : public FLAG {
// public:
//   BSSEEKERFLAG(const size_t f) : FLAG(f) {}
//   // pair-end:
//   //  if T-rich mate is +, then both mates are +;
//   //  if T-rich mate is -, then both mates are -;
//   // single-end:
//   //  0 for +; 16 for -.
//   bool is_Trich() const {
//     return FLAG::is_pairend() ? FLAG::is_Trich() : true;
//   }
//   bool is_Arich() const {
//     return FLAG::is_pairend() && FLAG::is_Arich();
//   }
//   bool is_revcomp() const {
//     if (FLAG::is_pairend())
//       return FLAG::is_revcomp() ? is_Trich() : is_Arich();
//     else
//       return FLAG::is_revcomp();
//   }
// };

// bool
// format_bsseeker_sam_record(sam_rec &samr) {

//   // string name, chrom, CIGAR, mate_name, seq, qual, orientation_str,
//   //   conversion_str, mismatch_str, mismatch_type_str, seq_genome_str;
//   // size_t flag, start, mapq_score, mate_start;
//   // int seg_len;

//   // std::istringstream iss(str);
//   // if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
//   //       >> mate_name >> mate_start >> seg_len >> seq >> qual
//   //       >> orientation_str >> conversion_str >> mismatch_str
//   //       >> mismatch_type_str >> seq_genome_str)) {
//   //   good = false;
//   //   throw runtime_error("malformed line in bs_seeker SAM format:\n" + str);
//   // }

//   // bs_seeker also doesn't keep sequencing quality information?
//   qual = string(seq.size(), 'h');

//   BSSEEKERFLAG Flag(flag);

//   samr.mr.r.set_name(name);
//   samr.mr.r.set_score(atoi(mismatch_str.substr(5).c_str()));
//   samr.mr.r.set_strand(Flag.is_revcomp() ? '-' : '+');

//   string new_seq, new_qual;
//   apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

//   if (Flag.is_revcomp()) {
//     revcomp_inplace(new_seq);
//     std::reverse(new_qual.begin(), new_qual.end());
//   }

//   samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size());
//   samr.mr.seq = new_seq;
//   samr.mr.qual = new_qual;

//   samr.is_Trich = Flag.is_Trich();
//   samr.is_mapping_paired = Flag.is_mapping_paired();

//   return good;
// }


// bool
// SAMReader::get_SAMRecord(const string &str, SAMRecord &samr) {
//   if (mapper == "bsmap")
//     return get_SAMRecord_bsmap(str, samr);
//   else if (mapper == "bismark")
//     return get_SAMRecord_bismark(str, samr);
//   else if (mapper == "bs_seeker")
//     return get_SAMRecord_bsseeker(str, samr);
//   else if (mapper == "general")
//     return get_SAMRecord_general(str, samr);
//   else
//     good = false;
//   return false;
// }


int
main(int argc, const char **argv) {

  try {

    string outfile;
    string mapper;
    size_t max_frag_len = 1000;
    size_t max_dangling = 500;
    size_t suffix_len = 1;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Convert the SAM/BAM output from "
                           "bismark or bs_seeker to MethPipe mapped read format",
                           "sam/bam_file");
    opt_parse.add_opt("output", 'o', "Name of output file",
                      false, outfile);
    opt_parse.add_opt("mapper", 'm',
                      "Original mapper: bismark, bs_seeker or general",
                      true, mapper);
    opt_parse.add_opt("suff", 's', "read name suffix length (default: 1)",
                      false, suffix_len);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (VERBOSE)
      cerr << "[input file: " << mapped_reads_file << "]" << endl
           << "[output file: "
           << (outfile.empty() ? "stdout" : outfile) << "]" << endl;

    // if (mapper == "bsmap")
    //   throw runtime_error("bsmap no longer supported [use general]");
    // if (mapper != "bismark" && mapper != "bs_seeker" && mapper != "general")
    //   throw runtime_error("mapper not supported: " + mapper);


    SAMReader sam_reader(mapped_reads_file, mapper);
    unordered_map<string, sam_rec> dangling_mates;

    size_t count = 0;
    sam_rec aln;
    while (sam_reader >> aln) {
      if (is_mapped(aln) && is_primary(aln)) {
        if (is_mapped_single_end(aln)) {
          if (is_a_rich(aln)) revcomp(aln);
          out << aln << '\n';
        }
        else { // is_mapped_paired(aln)
          const string read_name(remove_suffix(aln.qname, suffix_len));
          auto the_mate = dangling_mates.find(read_name);
          if (the_mate == end(dangling_mates)) {
            dangling_mates[read_name] = aln; // no mate yet
          }
          else { // found a mate
            if (is_t_rich(aln)) // earlier mate must have been t-rich
              swap(aln, the_mate->second);
            revcomp(aln); // put on same strand as earlier mate

            sam_rec merged;
            const size_t frag_len = merge_mates(suffix_len, max_frag_len,
                                                the_mate->second, aln, merged);
            if (frag_len <= max_frag_len)
              out << merged << '\n';
            else if (frag_len > 0)
              out << the_mate->second << '\n'
                  << aln << '\n';
            dangling_mates.erase(read_name);
          }

          // flush dangling_mates
          if (dangling_mates.size() > max_dangling) {
            unordered_map<string, sam_rec> to_keep;
            for (auto &&mates : dangling_mates)
              if (precedes_by_more_than(the_mate->second, aln, max_frag_len)) {
                if (is_a_rich(mates.second)) revcomp(mates.second);
                out << mates.second << endl;
              }
              else to_keep.insert(mates);
            swap(to_keep, dangling_mates);
          }
        }
        ++count;
      }
    }

    // flushing dangling_mates
    while (!dangling_mates.empty()) {
      if (is_a_rich(begin(dangling_mates)->second))
        revcomp(begin(dangling_mates)->second);
      out << begin(dangling_mates)->second << endl;
      dangling_mates.erase(begin(dangling_mates));
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
