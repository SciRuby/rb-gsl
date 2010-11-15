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
  @@mapping = [ 
    # numbered list
    #[/^(\s*)\((\d+)\)/, '\1\2.'],
    [/^(\s*)\((\d+)\)/, '1.'],
    # labelled list
    [/^(\s*):(\s*)(\S+)/, '\1[\2]'],
    # line contains a method. bulletize
    [/^(\s*)\-{3}\s/, '* \1'],
    # Emphasis
    [/\(\(\*/, '<b>'], [/\*\)\)/, '</b>'], # em?
    # Code
    [/\(\(\{/, '<tt>'], [/\}\)\)/, '</tt>'],
    # Variable
    [/\(\(\|/, '<tt>'], [/\|\)\)/, '</tt>'], # em?
    # Kbd input
    [/\(\(\%/, '<tt>'], [/\%\)\)/, '</tt>'],
    # Index, not supported.
    [/\(\(\:/, ''], [/\:\)\)/, ''],
    # Non-image hyperlink, auto detected by rdoc, use target=_top trick to open
    # external links in top level frame.
    [/\(\(\<(.*?)\|URL:http:(.*?)\>\)\)/, '{\1}[http:\2"target="_top]'],
    # Image hyperlink, auto detected by rdoc, use target=_top trick to open
    # external links in top level frame.
    [/\(\(\<(.*?)\|"IMG:http:(.*?)"\>\)\)/, '{\1}[http:\2]'],
    # Local hyperlink to top index
    # TODO: use link: and _top trick
    # Local HTML hyperlink, auto detected by rdoc
    # Add _rdoc suffix to stem
    [/\(\(\<(.*?)\|URL:(.*?).html([^>]*)\>\)\)/, '{\1}[link:files/rdoc/\2_rdoc.html\3]'],
    # Local non-HTML hyperlink, auto detected by rdoc
    [/\(\(\<(.*?)\|URL:(.*?)\>\)\)/, '{\1}[link:files/rdoc/\2]'],
    # Footnote, not supported.
    [/\(\(\-/, '(Note: '], [/\-\)\)/, ')'],
    # Verbatim text.
    [/\(\(\'/, '<tt>'], [/\'\)\)/, '</tt>'],
    # IRB prompt simplifiation
    [/irb\(.*?\):\d+:\d+[>*]/, '>>'],
  ]

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
    header_levels = [0]
    strip_mode = false
    method_block = false
    while line = @input.gets
      line.chomp!
      case @state
      when INSIDE_COMMENT
        if line =~ /^=end/
          line = ''
          @state=OUTSIDE_COMMENT
        else
          # Non-title heading lines
          if line =~ /^=(=+)\s+(.*)/
            header_levels += [0] until header_levels.length >= $1.length
            header_levels = header_levels[0, $1.length]
            header_levels[-1] += 1
            line = %Q{=#{$1} {}[link:index.html"name="#{header_levels.join('.')}] #{$2}}
          end
          # If method line, turn on "strip mode", which strips off 4 leading
          # spaces.
          if line =~ /^---\s/
            # Put rule line before blocks of methods
            @output.puts '# ---' if !method_block
            strip_mode = true
            method_block = true
          elsif line =~ /^=/
            # Put blank line after blocks of methods
            @output.puts '#' if method_block
            strip_mode = false
            method_block = false
          elsif strip_mode
            # Strip leading whitespace before mapping
            line.gsub!(/^\s{4}/,'  ')
            # Put blank line after blocks of methods
            @output.puts '#' if method_block
            method_block = false
          end
          @@mapping.each do |pattern,replacement|
            line.gsub!(pattern, replacement)
          end
        end
      when OUTSIDE_COMMENT
        if line =~ /^=begin/
          #line.gsub!(/^=begin/, '=begin rdoc')
          line = ''
          @state=INSIDE_COMMENT
        end
      else
        raise RdToRdocError.new("Unknown state #{@state}")
      end
      @output.print '#'
      @output.print ' ' unless line.empty?
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
