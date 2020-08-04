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
#include "bisulfite_utils.hpp"

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
using std::to_string;

// ADS: do we need to check that both mates are on the same strand? Or
// that they are on opposite strands?

static bool
is_mapped(const sam_rec &aln) {
  return !check_flag(aln, samflags::read_unmapped);
}

static bool
is_primary(const sam_rec &aln) {
  return !check_flag(aln, samflags::secondary_aln);
}

static bool
is_mapped_single_end(const sam_rec &aln) {
  return is_mapped(aln) &&
    (!check_flag(aln, samflags::read_paired) ||
     check_flag(aln, samflags::mate_unmapped));
}

inline bool
is_rc(const sam_rec &aln) {
  return check_flag(aln, samflags::read_rc);
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


static size_t
merge_mates(const size_t suffix_len, const size_t range,
            const sam_rec &one, const sam_rec &two, sam_rec &merged) {

  assert(is_rc(one) == false && is_rc(two) == true);
  if (one.pos > two.pos) {
    cerr << one << endl
         << two << endl
         << endl;
  }

  assert(is_t_rich(one) == is_t_rich(two));

  merged = one;

  const int one_s = one.pos;
  const int one_e = cigar_rseq_ops(one.cigar);

  const int two_s = two.pos;
  const int two_e = cigar_rseq_ops(two.cigar);

  const int spacer = two_s - one_e;
  if (spacer >= 0) {
    /* fragments longer enough that there is space between them: this
     * size of the spacer ("_") is determined based on the reference
     * positions of the two ends, and here we assume "one" maps to
     * positive genome strand.
     *
     * left                                                         right
     * one_s                    one_e      two_s                    two_e
     * [------------end1------------]______[------------end2------------]
     */
    merged.seq = one.seq + revcomp(two.seq);
    // ADS: need to take care of soft clipping in between;
    merged.cigar = one.cigar + to_string(spacer) + "N" + two.cigar;
    merged.qual = (one.qual == "*" ? one.qual : one.qual + revcomp(two.qual));

  }
  else {
    const int head = two_s - one_s;
    if (head >= 0) {
      /* fragment longer than or equal to the read length, but shorter
       * than twice the read length: this is determined by obtaining
       * the size of the "head" in the diagram below: the portion of
       * end1 that is not within [=]. If the read maps to the positive
       * strand, this depends on the reference start of end2 minus the
       * reference start of end1. For negative strand, this is
       * reference start of end1 minus reference start of end2.
       *
       * <======= head ================>
       *
       * left                                             right
       * one_s              two_s      one_e              two_e
       * [------------end1------[======]------end2------------]
       */
      truncate_cigar_r(merged.cigar, head);
      merged.cigar += two.cigar;
      const uint32_t merged_seq_len = cigar_qseq_ops(merged.cigar);
      // ADS: need to take care of soft clipping in between;
      merged.seq.resize(merged_seq_len);
      merged.seq += two.seq;
      if (merged.qual != "*") {
        merged.qual.resize(merged_seq_len);
        merged.qual += two.qual;
      }
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       * left                                       right
       * two_s            one_s    two_e            one_e
       * [--end2----------[============]----------end1--]
       */
      const int overlap = two_e - one_s;
      if (overlap >= 0) {
        truncate_cigar_r(merged.cigar, overlap);
        const uint32_t merged_seq_qlen = cigar_qseq_ops(merged.cigar);
        merged.seq.resize(merged_seq_qlen);
        if (merged.qual != "*")
          merged.qual.resize(merged_seq_qlen);
      }
    }
  }
  merged.rnext = "*";
  merged.pnext = 0;
  merged.tlen = 0;

  return two_e - one_s;
}

inline static bool
is_same_read(const size_t suffix_len, const sam_rec &a, const sam_rec &b) {
  return std::equal(begin(a.qname), end(a.qname) - suffix_len, begin(b.qname));
}

/********Above are functions for merging pair-end reads********/

static string
remove_suffix(const string &x, const size_t suffix_len) {
  return x.substr(0, x.size() - suffix_len);
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

// std::istringstream iss(str);
// if (!(iss
//       >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
//       >> mate_name >> mate_start >> seg_len >> seq >> qual
//       >> orientation_str >> conversion_str >> mismatch_str
//       >> mismatch_type_str >> seq_genome_str)) {
// good = false;
//   throw runtime_error("malformed line in bs_seeker SAM format:\n" + str);
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

    SAMReader sam_reader(mapped_reads_file);
    unordered_map<string, sam_rec> dangling_mates;

    out << sam_reader.get_header(); // includes newline

    size_t count = 0;
    sam_rec aln;

    while (sam_reader >> aln) {
      if (is_mapped(aln) && is_primary(aln)) {
        if (is_mapped_single_end(aln)) {
          if (is_a_rich(aln))
            revcomp(aln);
          out << aln << '\n';
        }
        else { // is_mapped_paired(aln)
          const string read_name(remove_suffix(aln.qname, suffix_len));
          auto the_mate = dangling_mates.find(read_name);
          if (the_mate == end(dangling_mates)) {
            dangling_mates[read_name] = aln; // no mate yet
          }
          else { // found a mate

            if (!is_rc(aln)) // earlier mate must have been a-rich
              swap(aln, the_mate->second);

            sam_rec merged;
            const int frag_len = merge_mates(suffix_len, max_frag_len,
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
                if (is_a_rich(mates.second))
                  revcomp(mates.second);
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
