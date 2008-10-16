#!/usr/bin/env ruby
# Convert RD documents to HTML

unless FileTest.directory?("html")
  puts("mkdir html/")
  Dir.mkdir("html")
end

dirs = ["rd"]
dirs.each do |dir|
  Dir.foreach(dir) do |f|
    if /\.rd$/ =~ f
      name = $`
      rd = "#{dir}/#{name}.rd"
      html = "html/#{name}.html"
      rdtime = File.mtime(rd)
      cmd = "rd2 -r rd/rd2html-ext-lib.rb --headline-secno --ref-extension --native-inline --head-element --with-part=head:head #{rd} > #{html}"
      if File.exist?(html)
        htmltime = File.mtime(html)
        a = rdtime <=> htmltime
        if a == 1
          puts cmd
          system "#{cmd}"
        end
      else
        puts cmd
        system "#{cmd}"
      end
    end
  end
end

