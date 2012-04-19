module Generator
  class Templater

    # Create a C file from a template. You must supply the name of the C file to write.
    #
    # If you wish to pass options to add by default, you may pass them in the second argument.
    # You may also supply a :header option, giving the filename of the output's header, if you so desire.
    def initialize output_filename, default_options = {}, &block
      @output_filename = output_filename
      @header_filename = default_options.delete(:header)
      @default_options = default_options
      @declarations    = []
      @begun           = false
      yield self
    end

    def header name, options = {}
      raise("please provide header before templates") if @begun
      add_internal(name, false, options)
    end

    # Add a .template.c file to the current output. Takes templates as mandatory argument (e.g., for numbmm.template.c
    # and transp.template.c, give it %w{numbmm transp}).
    #
    # Options
    # =======
    #  * :in  :: subdirectory to search in
    def template templates, options = {}
      # Header has to go first
      header(@header_filename, options) unless @header_filename.nil?
      @begun = true

      templates = [templates] if templates.is_a?(String)

      @declarations = @declarations.concat(add_internal(templates, true, options))
    end

  protected
    def add_internal names, as_template, options
      opts = options.dup
      dir = opts.delete(:in) || @default_options[:in]

      if as_template
        template_hash = @default_options.merge(opts).select { |k,v| k.to_s =~ /^[A-Z_0-9]+$/}
        apply_each_template_type(@output_filename, names, dir, template_hash)
      else
        # Pretend names is just one name for this part:
        template_path = File.join(relative_ext_path(dir), "#{names}.template.c")
        `cat #{template_path} > ./#{@output_filename}`
        []
      end
    end

    def relative_ext_path subdir=nil
      subdir.nil? || subdir.size == 0 ? File.join($RELATIVE_PATH, "ext", "nmatrix") : File.join($RELATIVE_PATH, "ext", "nmatrix", subdir)
    end

    # Returns an array of C function declarations, which can then be put in nmatrix.h
    def apply_each_template_type output_path, names, dir, type_hash
      raise(NotImplementedError, "Too many keys, please fix") if type_hash.keys.size > 2

      template_path = !dir.nil? && dir.size > 0 ? File.join(relative_ext_path, dir) : File.join($RELATIVE_PATH, "ext", "nmatrix")

      declarations = []

      a = type_hash.keys[0]
      b = type_hash.keys[1] if type_hash.size > 1

      raise(ArgumentError, "no template keys found") if !defined?(a) || a.nil?

      type_hash[a].each do |type_a|

        # Currently allows just one or two types. TODO: Fix this so it allows as many types as we could want. Low priority.
        if defined?(b) && !b.nil?
          type_hash[b].each do |type_b|
            names.each do |name|
              full_template_path = File.join(template_path, "#{name}.template.c")
              raise(ArgumentError, "template not found at #{full_template_path} relative to #{Dir.pwd}") unless File.exists?(full_template_path)

              declarations << Generator.template(full_template_path, output_path, {a => type_a, b => type_b})
            end
          end
        else
          names.each do |name|

            full_template_path = File.join(template_path, "#{name}.template.c")
            raise(ArgumentError, "template not found at #{full_template_path} relative to #{Dir.pwd}") unless File.exists?(full_template_path)

            declarations << Generator.template(full_template_path, output_path, {a => type_a})
          end
        end

      end

      declarations
    end
  end
end