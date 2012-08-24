#!/usr/bin/ruby

# A helper file for generating and maintaining template tables.

def nullify(disabled)
	DTYPES.map { |t| if disabled.include?(t) then :NULL else t end }
end

DTYPES = [
	:uint8_t,
	:int8_t,
	:int16_t,
	:int32_t,
	:int64_t,
	:float32_t,
	:float64_t,
	:Complex64,
	:Complex128,
	:Rational32,
	:Rational64,
	:Rational128,
	:RubyObject
]

EWOPS = [
	:EW_ADD,
	:EW_SUB,
	:EW_MUL,
	:EW_DIV,
	:EW_MOD
]

LR_ALLOWED = {
	:uint8_t	 		=> nullify([:RubyObject]),
	:int8_t				=> nullify([:RubyObject]),
	:int16_t			=> nullify([:RubyObject]),
	:int32_t			=> nullify([:RubyObject]),
	:int64_t			=> nullify([:RubyObject]),
	:float32_t		=> nullify([:RubyObject]),
	:float64_t		=> nullify([:RubyObject]),
	:Complex64		=> nullify([:RubyObject]),
	:Complex128		=> nullify([:RubyObject]),
	:Rational32		=> nullify([:float32_t, :float64_t, :Complex64, :Complex128, :RubyObject]),
	:Rational64		=> nullify([:float32_t, :float64_t, :Complex64, :Complex128, :RubyObject]),
	:Rational128	=> nullify([:float32_t, :float64_t, :Complex64, :Complex128, :RubyObject]),
	:RubyObject		=> nullify(DTYPES - [:RubyObject])
}

lines =
case ARGV[0]
when 'OPLR'
	'{' +
	EWOPS.map do |op|
	
		'{' +
		DTYPES.map do |l_dtype|
		
			'{' +
			LR_ALLOWED[l_dtype].map do |r_dtype|
				if r_dtype == :NULL
					'NULL'
				else
					"fun<#{op}, #{l_dtype}, #{r_dtype}>"
				end
			end.join(', ') +
			'}'
		
		end.join(",\n") +
		'}'
	
	end.join(",\n") +
	'}'

when 'LR'
	'{' +
		DTYPES.map do |l_dtype|
		
			'{' +
			LR_ALLOWED[l_dtype].map do |r_dtype|
				if r_dtype == :NULL
					'NULL'
				else
					"fun<#{l_dtype}, #{r_dtype}>"
				end
			end.join(', ') +
			'}'
		
		end.join(",\n") +
		'}'
end

puts lines
