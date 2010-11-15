#!/usr/local/bin/ruby -w
#
# A quick and dirty converter between RD (rdtool) and RDoc formats.
#
# Author:: Hugh Sasee
# Created:: 19-DEC-2005
# Last Modified:: 19-DEC-2005
#
#

class RdToRdocError < RuntimeError; end

class RdToRdoc
  OUTSIDE_COMMENT = :OUTSIDE_COMMENT
  INSIDE_COMMENT = :INSIDE_COMMENT
  @@mapping = { 
    /^(\s*)\((\d+)\)/ => '\1\2', # numbered list
    /^(\s*):(\s*)(\S+)/ => '\1[\2]', # labelled list
    /^(\s*)\-{3}/ => '\1   ', # line contains a method. Rdoc should cope
    /\(\(\*(.*?)\*\)\)/ => '<em>\1</em>',
    /\(\(\{(.*?)\}\)\)/ => '<tt>\1</tt>',
    /\(\(\|(.*?)\|\)\)/ => '#\1',
    /\(\(\%(.*?)\%\)\)/ => '<tt>\1</tt>', # Kbd input
    /\(\(\:(.*?)\:\)\)/ => '\1', # Index, not supported.
    /\(\(\<(.*?)\>\)\)/ => '\1', # Hyperlink, auto detected by rdoc
    /\(\(\-(.*?)\-\)\)/ => '(Note: \1)', # Footnote, not supported.
    /\(\(\'(.*?)\'\)\)/ => '(Note: \1)'  # Verbatim text.
  }

  # The constructor takes an input object that responds to readline,
  # and an output object that responds to puts.  Use is made of
  # chomp when reading to allow conversion between line-endings if
  # necessary.  These are not specified here, so it is up to the
  # user to set $/ and $\ correctly.  IO-like objects are used for
  # greatest flexibility,
  def initialize(in_io, out_io)
    @state = OUTSIDE_COMMENT
    @verbose = false
    @input = in_io
    @output = out_io
  end

  def process()
    while line = @input.gets
      case @state
      when INSIDE_COMMENT
        if line =~ /^=end/
          @state=OUTSIDE_COMMENT
        else
          @@mapping.each do |key,value|
            line.gsub!(key, value)
          end
        end
      when OUTSIDE_COMMENT
        if line =~ /^=begin/
          line.gsub!(/^=begin/, '=begin rdoc')
          @state=INSIDE_COMMENT
        end
      else
        raise RdToRdocError.new("Unknown state #{@state}")
      end
      @output.puts line
    end
    @input.close
    @output.close
  end
end

if __FILE__ == $0
  case ARGV.size
  when 0
    input = $stdin
    output = $stdout
  when 1
    input = open(ARGV[0]) or raise "unable to open #{ARGV[0]} for reading #{$!}"
    output = $stdout
  when 2
    input = open(ARGV[0]) or raise "unable to open #{ARGV[0]} for reading #{$!}"
    output = open(ARGV[1], "w") or raise "unable to open #{ARGV[0]} for writing #{$!}"
  else
    raise "usage #{$0} [input] [output]"
  end

  RdToRdoc.new(input, output).process()
end
