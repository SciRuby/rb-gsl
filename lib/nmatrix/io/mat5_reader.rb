# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == io/matlab/mat5_reader.rb
#
# Matlab version 5 .mat file reader (and eventually writer too).
#


require_relative 'mat_reader.rb'

module NMatrix::IO::Matlab
	# Reader (and eventual writer) for a version 5 .mat file.
	class Mat5Reader < MatReader
		attr_reader :file_header, :first_tag_field, :first_data_field

    class Compressed
      include Packable
      # include TaggedDataEnumerable

      attr_reader :byte_order

      def initialize(stream = nil, byte_order = nil, content_or_bytes = nil)
        @stream			= stream
        @byte_order	= byte_order

        if content_or_bytes.is_a?(String)
          @content = content_or_bytes

        elsif content_or_bytes.is_a?(Fixnum)
          @padded_bytes = content_or_bytes
        #else
        #  raise ArgumentError, "Need a content string or a number of bytes; content_or_bytes is #{content_or_bytes.class.to_s}."
        end
      end

      def compressed
        require "zlib"
        # [2..-5] removes headers
        @compressed ||= Zlib::Deflate.deflate(content)
      end

      def content
        @content ||= extract
      end

      def padded_bytes
        @padded_bytes ||= content.size % 4 == 0 ? content.size : (content.size / 4 + 1) * 4
      end

      def write_packed(packedio, options)
        packedio << [compressed, {:bytes => padded_bytes}.merge(options)]
      end

      def read_packed(packedio, options)
        @compressed = (packedio >> [String, options]).first
        content
      end

      protected
      def extract
        require 'zlib'

        zstream = Zlib::Inflate.new #(-Zlib::MAX_WBITS) # No header

        returning(zstream.inflate(@compressed)) do
          zstream.finish
          zstream.close
        end
      end
    end

    MatrixDataStruct = Struct.new(
      :cells, :logical, :global, :complex, :nonzero_max,
      :matlab_class, :dimensions, :matlab_name, :real_part,
      :imaginary_part, :row_index, :column_index)

    class MatrixData < MatrixDataStruct
      include Packable

      def write_packed(packedio, options)
        raise NotImplementedError
        packedio << [info, {:bytes => padded_bytes}.merge(options)]
      end

      # Figure out the appropriate Ruby type to convert to, and do it. There are basically two possible types: NMatrix
      # and Ruby Array. This function is recursive, so an Array is going to contain other Arrays and/or NMatrix objects.
      #
      # mxCELL types (cells) will be converted to the Array type.
      #
      # mxSPARSE and other types will be converted to NMatrix, with the appropriate stype (:yale or :dense, respectively).
      #
      # See also to_nm, which is responsible for NMatrix instantiation.
      def to_ruby
        case matlab_class
        when :mxSPARSE	then	return to_nm
        when :mxCELL		then	return self.cells.collect { |c| c.to_ruby }
        else									return to_nm
        end
      end

      # Try to determine what dtype and such to use.
      #
      # TODO: Needs to be verified that unsigned MATLAB types are being converted to the correct NMatrix signed dtypes.
      def guess_dtype_from_mdtype
        dtype = MatReader::MDTYPE_TO_DTYPE[self.real_part.tag.data_type]

        return dtype unless self.complex

        dtype == :float32 ? :complex64 : :complex128
      end

      # Unpacks and repacks data into the appropriate format for NMatrix.
      #
      # If data is already in the appropriate format, does not unpack or repack, just returns directly.
      #
      # Complex is always unpacked and repacked, as the real and imaginary components must be merged together (MATLAB
      # stores them separately for some crazy reason).
      def repacked_data(to_dtype = nil)
        real_mdtype = self.real_part.tag.data_type

        if !self.complex and MatReader::NO_REPACK.include?(real_mdtype)
          self.real_part.data

        else
          to_dtype ||=
          begin
            imag_mdtype = self.imaginary_part.tag.data_type

            # Figure out what dtype we need to convert to and get arguments for pack
            if self.complex && real_mdtype != imag_mdtype
              # TODO: Expose upcasting in NMatrix namespace and allow this to take different dtypes in real and imag components.
              raise NotImplementedError, 'MATLAB dtypes differ for real and imaginary components, please specify dtype as argument'
            end

            MatReader::MDTYPE_TO_DTYPE[real_mdtype]
          end

          repack_args = MatReader::DTYPE_PACK_ARGS[to_dtype]

          unpacked_data(real_mdtype, imag_mdtype).pack(repack_args)
        end
      end

      def unpacked_data(real_mdtype = nil, imag_mdtype = nil)
        # Get Matlab data type and unpack args
        real_mdtype			||= self.real_part.tag.data_type
        real_unpack_args	= MatReader::MDTYPE_UNPACK_ARGS[real_mdtype]

        # zip real and complex components together, or just return real component
        if self.complex
          imag_mdtype			||= self.imaginary_part.tag.data_type
          imag_unpack_args	= MatReader::MDTYPE_UNPACK_ARGS[imag_mdtype]

          unpacked_real = self.real_part.data.unpack(real_unpack_args)
          unpacked_imag = self.imaginary_part.data.unpack(imag_unpack_args)

          unpacked_real.zip(unpacked_imag).flatten

        else
          self.real_part.data.unpack(real_unpack_args)
        end
      end


      # Unpacks and repacks index data into the appropriate format for NMatrix.
      #
      # If data is already in the appropriate format, does not unpack or repack, just returns directly.
      def repacked_indices(to_itype)
        return [row_index.data, column_index.data] if to_itype == :uint32 # No need to re-pack -- already correct

        unpack_args = [MatReader::MDTYPE_UNPACK_ARGS[:miINT32]] * (self.row_index.data.length / 4)
        unpacked_row_indices = self.row_index.data.unpack(*unpack_args)
        repacked_row_indices = unpacked_row_indices.pack(*([MatReader::ITYPE_PACK_ARGS[to_itype][1]] * unpacked_row_indices.size))

        unpack_args = [MatReader::MDTYPE_UNPACK_ARGS[:miINT32]] * (self.column_index.data.length / 4)
        unpacked_column_indices = self.column_index.data.unpack(*unpack_args)
        repacked_column_indices = unpacked_column_indices.pack(*([MatReader::ITYPE_PACK_ARGS[to_itype][1]] * unpacked_column_indices.size))

        [repacked_row_indices, repacked_column_indices]
      end


      # Create an NMatrix from a MATLAB .mat (v5) matrix.
      #
      # This function matches the storage type exactly. That is, a regular matrix in MATLAB will be a dense NMatrix, and
      # a sparse (old Yale) matrix in MATLAB will be a :yale (new Yale) matrix in NMatrix.
      #
      # Note that NMatrix has no old Yale type, so this uses a semi-hidden version of the NMatrix constructor to pass in
      # -- as directly as possible -- the stored bytes in a MATLAB sparse matrix. This constructor should also be used
      # for other IO formats that want to create sparse matrices from IA and JA vectors (e.g., SciPy).
      #
      # This is probably not the fastest code. An ideal solution would be a C plugin of some sort for reading the MATLAB
      # .mat file. However, .mat v5 is a really complicated format, and lends itself to an object-oriented solution.
      def to_nm(dtype = nil)
        # Hardest part is figuring out from_dtype, from_index_dtype, and dtype.
        dtype			||= guess_dtype_from_mdtype
        from_dtype	= MatReader::MDTYPE_TO_DTYPE[self.real_part.tag.data_type]

        # Create the same kind of matrix that MATLAB saved.
        case matlab_class
        when :mxSPARSE
          raise(NotImplementedError, "expected .mat row indices to be of type :miINT32") unless row_index.tag.data_type == :miINT32
          raise(NotImplementedError, "expected .mat column indices to be of type :miINT32") unless column_index.tag.data_type == :miINT32

          to_itype = NMatrix.itype_by_shape(dimensions)

          #require 'pry'
          #binding.pry

          # MATLAB always uses :miINT32 for indices according to the spec
          NMatrix.new(:yale, dimensions, dtype, *repacked_indices(to_itype), repacked_data(dtype), from_dtype)

        else
          # Call regular dense constructor.
          NMatrix.new(:dense, dimensions, unpacked_data, dtype)
        end
      end

      def read_packed(packedio, options)
        flags_class, self.nonzero_max = packedio.read([Element, options]).data

        self.matlab_class   = MatReader::MCLASSES[flags_class % 16]
        #STDERR.puts "Matrix class: #{self.matlab_class}"

        self.logical        = (flags_class >> 8) % 2 == 1 ? true : false
        self.global         = (flags_class >> 9) % 2 == 1 ? true : false
        self.complex        = (flags_class >> 10) % 2 == 1 ? true : false
        #STDERR.puts "nzmax: #{self.nonzero_max}"

        dimensions_tag_data = packedio.read([Element, options])
        self.dimensions     = dimensions_tag_data.data
        #STDERR.puts "dimensions: #{self.dimensions}"

        begin
          name_tag_data			= packedio.read([Element, options])
          self.matlab_name	= name_tag_data.data.is_a?(Array) ? name_tag_data.data.collect { |i| i.chr }.join('') : name_tag_data.data.chr

        rescue ElementDataIOError => e
          STDERR.puts "ERROR: Failure while trying to read Matlab variable name: #{name_tag_data.inspect}"
          STDERR.puts 'Element Tag:'
          STDERR.puts "    #{e.tag}"
          STDERR.puts 'Previously, I read these dimensions:'
          STDERR.puts "    #{dimensions_tag_data.inspect}"
          STDERR.puts "Unpack options were: #{options.inspect}"
          raise(e)
        end

        #STDERR.puts [flags_class.to_s(2), self.complex, self.global, self.logical, nil, self.mclass, self.nonzero_max].join("\t")
        if self.matlab_class == :mxCELL
          # Read what may be a series of matrices
          self.cells = []
          STDERR.puts("Warning: Cell array does not yet support reading multiple dimensions") if dimensions.size > 2 || (dimensions[0] > 1 && dimensions[1] > 1)
          number_of_cells = dimensions.inject(1) { |prod,i| prod * i }
          number_of_cells.times { self.cells << packedio.read([Element, options]) }

        else
          read_opts = [RawElement, {:bytes => options[:bytes], :endian => :native}]

          if self.matlab_class == :mxSPARSE
            self.column_index = packedio.read(read_opts)
            self.row_index    = packedio.read(read_opts)

            # STDERR.puts "row and col indices: #{self.row_index.inspect}, #{self.column_index.inspect}"
          end

          self.real_part			= packedio.read(read_opts)
          self.imaginary_part	= packedio.read(read_opts) if self.complex
        end
      end

      def ignore_padding packedio, bytes
        packedio.read([Integer, {:unsigned => true, :bytes => bytes}]) if bytes > 0
      end
    end


		MDTYPE_UNPACK_ARGS =
		MatReader::MDTYPE_UNPACK_ARGS.merge({
			:miCOMPRESSED	=> [Compressed, {}],
			:miMATRIX			=> [MatrixData, {}]
		})
		# include TaggedDataEnumerable
		
		FIRST_TAG_FIELD_POS = 128
		
		####################
		# Instance Methods #
		####################

		def initialize(stream, options = {})
			super(stream, options)
			@file_header = seek_and_read_file_header
		end

		def to_a
			returning(Array.new) do |ary|
				self.each { |el| ary << el }
			end
		end

		def to_ruby
			ary = self.to_a
			
			if ary.size == 1
				ary.first.to_ruby 
			else
				ary.collect { |item| item.to_ruby }
			end
		end

		def guess_byte_order
			stream.seek(Header::BYTE_ORDER_POS)
			mi = stream.read(Header::BYTE_ORDER_LENGTH)
			stream.seek(0)
			mi == 'IM' ? :little : :big
		end

		def seek_and_read_file_header
			stream.seek(0)
			stream.read(FIRST_TAG_FIELD_POS).unpack(Header, {:endian => byte_order})
		end

		def each(&block)
			stream.each(Element, {:endian => byte_order}) do |element|
				if element.data.is_a?(Compressed)
					StringIO.new(element.data.content, 'rb').each(Element, {:endian => byte_order}) do |compressed_element|
						yield compressed_element.data
					end
					
				else
					yield element.data
				end
			end
			
			# Go back to the beginning in case we want to do it again.
			stream.seek(FIRST_TAG_FIELD_POS)
			
			self
		end
		
		####################
		# Internal Classes #
		####################
		
		class Header < Struct.new(:desc, :data_offset, :version, :endian)

			include Packable

			BYTE_ORDER_LENGTH		= 2
			DESC_LENGTH					= 116
			DATA_OFFSET_LENGTH	= 8
			VERSION_LENGTH			= 2
			BYTE_ORDER_POS			= 126

			## TODO: TEST WRITE.
			def write_packed(packedio, options)
				packedio <<	[desc,				{:bytes => DESC_LENGTH				}] <<
                    [data_offset,	{:bytes => DATA_OFFSET_LENGTH	}] <<
                    [version,			{:bytes => VERSION_LENGTH			}] <<
                    [byte_order,	{:bytes => BYTE_ORDER_LENGTH	}]
			end

			def read_packed(packedio, options)
				self.desc, self.data_offset, self.version, self.endian = packedio >>
          [String,	{:bytes => DESC_LENGTH																	}] >>
					[String,	{:bytes => DATA_OFFSET_LENGTH														}] >>
					[Integer,	{:bytes => VERSION_LENGTH, :endian => options[:endian]	}] >>
					[String,	{:bytes => 2																						}]
				
				self.desc.strip!
				self.data_offset.strip!
				self.data_offset = nil if self.data_offset.empty?
				
				self.endian == 'IM' ? :little : :big
			end
		end
		
		class Tag < Struct.new(:data_type, :raw_data_type, :bytes, :small)
			include Packable
			
			DATA_TYPE_OPTS = BYTES_OPTS = {:bytes => 4, :signed => false}
			LENGTH = DATA_TYPE_OPTS[:bytes] + BYTES_OPTS[:bytes]
			
			## TODO: TEST WRITE.
			def write_packed packedio, options
				packedio << [data_type, DATA_TYPE_OPTS] << [bytes, BYTES_OPTS]
			end
			
			def small?
				self.bytes > 0 and self.bytes <= 4
			end
			
			def size
				small? ? 4 : 8
			end
			
			def read_packed packedio, options
				self.raw_data_type = packedio.read([Integer, DATA_TYPE_OPTS.merge(options)])
				
				# Borrowed from a SciPy patch
				upper = self.raw_data_type >> 16
				lower = self.raw_data_type & 0xFFFF
				
				if upper > 0
					# Small data element format
					raise IOError, 'Small data element format indicated, but length is more than 4 bytes!' if upper > 4
					
					self.bytes					= upper
					self.raw_data_type	= lower
					
				else
					self.bytes = packedio.read([Integer, BYTES_OPTS.merge(options)])
				end
				
				self.data_type = MatReader::MDTYPES[self.raw_data_type]
			end

			def inspect
				"#<#{self.class.to_s} data_type=#{data_type}[#{raw_data_type}][#{raw_data_type.to_s(2)}] bytes=#{bytes} size=#{size}#{small? ? ' small' : ''}>"
			end
		end


		class ElementDataIOError < IOError
			attr_reader :tag
			
			def initialize(tag = nil, msg = nil)
				@tag = tag
				super msg
			end

			def to_s
				@tag.inspect + "\n" + super
			end
		end


		class Element < Struct.new(:tag, :data)
			include Packable

			def write_packed packedio, options
				packedio << [tag, {}] << [data, {}]
			end

			def read_packed(packedio, options)
				raise(ArgumentError, 'Missing mandatory option :endian.') unless options.has_key?(:endian)
				
				tag = packedio.read([Tag, {:endian => options[:endian]}])
				#STDERR.puts tag.inspect
				data_type = MDTYPE_UNPACK_ARGS[tag.data_type]
				
				self.tag = tag
				#STDERR.puts self.tag.inspect
				
				raise ElementDataIOError.new(tag, "Unrecognized Matlab type #{tag.raw_data_type}") if data_type.nil?

				if tag.bytes == 0
					self.data = []
					
				else
					number_of_reads = data_type[1].has_key?(:bytes) ? tag.bytes / data_type[1][:bytes] : 1
					data_type[1].merge!({:endian => options[:endian]})

					if number_of_reads == 1
						self.data = packedio.read(data_type)
						
					else
						self.data =
						returning(Array.new) do |ary|
							number_of_reads.times { ary << packedio.read(data_type) }
						end
					end
				
					begin
						ignore_padding(packedio, (tag.bytes + tag.size) % 8) unless [:miMATRIX, :miCOMPRESSED].include?(tag.data_type)
						
					rescue EOFError
						STDERR.puts self.tag.inspect
						raise(ElementDataIOError.new(tag, "Ignored too much"))
					end
				end
			end

			def ignore_padding(packedio, bytes)
				if bytes > 0
					#STDERR.puts "Ignored #{8 - bytes} on #{self.tag.data_type}"
					ignored = packedio.read(8 - bytes)
					ignored_unpacked = ignored.unpack("C*")
					raise(IOError, "Nonzero padding detected: #{ignored_unpacked}") if ignored_unpacked.any? { |i| i != 0 }
				end
			end

			def to_ruby
				data.to_ruby
			end
		end

		# Doesn't unpack the contents of the element, e.g., if we want to handle
		# manually, or pass the raw string of bytes into NMatrix.
		class RawElement < Element
			def read_packed(packedio, options)
				raise(ArgumentError, 'Missing mandatory option :endian.') unless options.has_key?(:endian)

				self.tag	= packedio.read([Tag,			{:endian => options[:endian]											}])
				self.data	= packedio.read([String,	{:endian => options[:endian], :bytes => tag.bytes	}])

				begin
					ignore_padding(packedio, (tag.bytes + tag.size) % 8) unless [:miMATRIX, :miCOMPRESSED].include?(tag.data_type)
					
				rescue EOFError
					STDERR.puts self.tag.inspect
					raise ElementDataIOError.new(tag, 'Ignored too much.')
				end
			end
		end
		
		#####################
		# End of Mat5Reader #
		#####################
		
	end
end
