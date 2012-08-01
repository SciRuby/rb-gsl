#!/usr/local/bin/ruby
# Generate LIR_DTYPE_TEMPLATE_TABLE macro.
ary = %w{uint8_t int8_t int16_t int32_t int64_t float32_t float64_t Complex64 Complex128 Rational32 Rational64 Rational128 RubyObject}
bry = ary.clone
iry = %w{uint8_t uint16_t uint32_t uint64_t}

lines = []
lines << "#define LRI_DTYPE_TEMPLATE_TABLE(fun, ret, ...)"
lines << "static ret (*ttable[NUM_DTYPES][NUM_DTYPES][NUM_ITYPES])(__VA_ARGS__) = {"

curr_line = ""

ary.each do |a|
  lines << curr_line if curr_line.size > 0
  curr_line = "{"

  bry.each do |b|
    curr_line += "{"

    curr_line_list = []
    iry.each do |i|
      curr_line_list << "fun<#{a},#{b},#{i}>"
    end

    curr_line += curr_line_list.join(',') + "},"
  end

  curr_line += "},"

  lines << curr_line
end

puts lines.join("\\\n") + "};"
