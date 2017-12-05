#!/usr/bin/env ruby
# ucsc2vcf.rb
# programmed by (c) MISHIMA, Hiroyuki (hmishima at nagasaki-u.ac.jp)
# License: the MIT/X11 license

#Version = "20120202"
Version = "20171206"

require 'optparse'
require 'bio'
require 'bio-ucsc-api'
require 'striuct'

module Ucsc2Vcf
  class Columns < Striuct.new
    member :bin, Integer, &->v{Integer v}
    member :chrom, String
    member :chromStart, Integer, &->v{Integer v}
    member :chromEnd, Integer, &->v{Integer v}
    member :name, String
    member :score, Integer, &->v{Integer v}
    member :strand, String
    member :refNCBI, String
    member :refUCSC, String
    member :observed, String
    member :molType, String
    member :klass, String
    member :valid, String
    member :avHet, Float, &->v{Float v}
    member :avHetSE, Float, &->v{Float v}
    member :func, String
    member :locType, String
    member :weight, Integer, &->v{Integer v}
    member :exceptions, String
    member :submitterCount, Integer, &->v{Integer v}
    member :submitters, String
    member :alleleFreqCount, Integer, &->v{Integer v}
    member :alleles, String
    member :alleleNs, String
    member :alleleFreqs, String
    member :bitfields, String
  end

  class VCFrecord < Striuct.new
    member  :chrom, String
    member  :pos, Integer
    member  :id, String
    member  :ref, String
    member  :alt, String
    member  :qual, String
    default :qual, "."
    member  :filter, String
    default :filter, "."
    member  :info, String
    member  :col_format, String
    default :col_format, "GT"
    member  :samplegt, String
    default :samplegt, "0/1"
    
    def to_s
      [ chrom, pos, id, ref, alt, qual,
        filter, info, col_format, samplegt].join("\t")
    end
  end

  class Application
    def initialize(opts)
      @opts = opts
      if opts[:indel]
        load_reference
      end
    end

    def output_header
      puts <<"EOF"
##fileformat=VCFv4.1
##source=UCSC_snp135
##dbSNP_BUILD_ID=135
##reference=GRCh37_hg19
##program=ucsc2vcf.rb
##INFO=<ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
##INFO=<ID=CLASS,Number=1,Type=String,Description="UCSC variation class">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdbSNP135
EOF
    end

    def load_reference
      @reference = Bio::Ucsc::Reference.load(@opts[:ref])
    end

    def refseq(pos)
      @reference.find_by_interval(Bio::GenomicInterval.parse(pos))
    end

    def locus_single(cols)
      vcf = VCFrecord.new
      vcf.chrom = cols.chrom
      vcf.pos = cols.chromStart + 1
      vcf.id = cols.name
      case cols.strand
      when "+"
        vcf.ref = cols.refUCSC
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==cols.refUCSC}.join(",")
      when "-"
        vcf.ref = cols.refUCSC
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==revcmp(cols.refUCSC)} \
          .map{|x|revcmp(x)}.join(",")
      else
        raise "Should not happen!"
      end
      vcf.info = "RSPOS=#{Integer(cols.chromStart)+1};CLASS=#{cols.klass}"
      puts vcf.to_s
    end

    def locus_insertion(cols)
      vcf = VCFrecord.new
      vcf.chrom = cols.chrom
      vcf.pos = cols.chromStart
      vcf.id = cols.name
      previous = refseq("#{cols.chrom}:#{Integer(cols.chromStart)}")
      vcf.ref = previous
      case cols.strand
      when "+"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==cols.refUCSC} \
          .map{|x|if x == "-" then previous else previous + x end}.join(",")
      when "-"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==revcmp(cols.refUCSC)} \
          .map{|x|revcmp(x)} \
          .map{|x|if x == "-" then previous else previous + x end}.join(",")
      else
        raise "Should not happen!"
      end
      vcf.info = "RSPOS=#{Integer(cols.chromStart)+1};CLASS=#{cols.klass}"
      puts vcf
    end

    def locus_deletion(cols)
      vcf = VCFrecord.new
      vcf.chrom = cols.chrom
      vcf.pos = cols.chromStart
      vcf.id = cols.name
      previous = refseq("#{cols.chrom}:#{Integer(cols.chromStart)}")
      vcf.ref = previous + cols.refUCSC
      case cols.strand
      when "+"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==cols.refUCSC} \
          .map{|x|if x == "-" then previous else previous + x end}.join(",")
      when "-"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==revcmp(cols.refUCSC)} \
          .map{|x|revcmp(x)} \
          .map{|x|if x == "-" then previous else previous + x end}.join(",")
      else
        raise "Should not happen!"
      end
      vcf.info = "RSPOS=#{Integer(cols.chromStart)+1};CLASS=#{cols.klass}"
      puts vcf
    end

    def locus_substitution(cols)
      vcf = VCFrecord.new
      vcf.chrom = cols.chrom
      vcf.pos = Integer(cols.chromStart) + 1
      vcf.id = cols.name
      vcf.ref = cols.refUCSC
      case cols.strand
      when "+"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==cols.refUCSC}.join(",")
      when "-"
        vcf.alt =
          cols.observed.split("/") \
          .reject{|x|x==revcmp(cols.refUCSC)} \
          .map{|x|revcmp(x)}.join(",")
      else
        raise "Should not happen!"
      end
      vcf.info = "RSPOS=#{Integer(cols.chromStart)+1};CLASS=#{cols.klass}"
      puts vcf
    end

    def locus(cols)
      case
      when (/\A[ACGT]\z/ =~ cols.refUCSC &&
            cols.observed.split("/").all?{|x|/\A[ACGT]\z/ =~ x})
        locus_single(cols) if @opts[:snv]
      when (cols.refUCSC == "-" &&
            cols.observed.split("/").none?{|x|/[^ACGT-]/ =~ x})
        locus_insertion(cols) if @opts[:indel]
      when (/[^ACGT]/ !~ cols.refUCSC &&
            cols.observed.split("/").none?{|x|/[^ACGT-]/ =~ x} &&
            cols.observed.split("/").any?{|x|x == "-"})
        locus_deletion(cols) if @opts[:indel]
      when (/[^ACGT]/ !~ cols.refUCSC &&
            cols.observed.split("/").none?{|x|/[^ACGT-]/ =~ x})
        locus_substitution(cols) if @opts[:indel]
      else
        # large indels are ignored
        # $stderr.puts cols.inspect
        nil
      end
    end

    def revcmp(na)
      Bio::Sequence::NA.new(na).reverse_complement.upcase
    end

    def run
      output_header
      ARGF.each_line.each do |row|
        locus Columns.new(*(row.chomp.split("\t")))
      end
    end
  end
end

if __FILE__ == $0
  opts = Hash.new
  opts[:samples] = Hash.new
  
  # ARGV[0] = "--help" if ARGV.length == 0
  ARGV.options do |o|
    o.banner = "ucsc2vcf.rb [option] [filename(or stdin)]"
    o.on('-s', '--snv',
         'output SNV/SNP loci (can be combined with -i') do
      opts[:snv] = true
    end      
    o.on('-i', '--indel',
         'output SNV/SNP loci (can be combined with -s') do
      opts[:indel] = true
    end
    o.on('-r filename', '--reference', String,
         'reference in the 2bit file format (required unless -s option is added)') do |x|
      opts[:ref] = x
    end
    o.separator "    -v, --version                    show version information"
    o.separator "    -h, --help                       show this message"
    o.separator " last update: #{o.version}"
    o.parse!
  end  
  Ucsc2Vcf::Application.new(opts).run
end
