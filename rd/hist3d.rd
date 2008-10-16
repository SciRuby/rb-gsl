=begin
= GSL::Histogram3d class

== Class methods
--- GSL::Histogram3d.alloc(nx, ny, nz)
--- GSL::Histogram3d.alloc(xrange, yrange, zrange)
--- GSL::Histogram3d.alloc(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax)
--- GSL::Histogram3d.alloc(nx, [xmin, xmax], ny, [ymin, ymax], nz, [zmin, zmax])
    Constructors

--- GSL::Histogram3d.memcpy(h1, h2)
--- GSL::Histogram3d.equal_bins_p(h1, h2)
--- GSL::Histogram3d.equal_bins_?p(h1, h2)

== Methods
--- GSL::Histogram3d#set_ranges(xrange, yrange, zrange)
--- GSL::Histogram3d#set_ranges_uniform(xmin, xmax, ymin, ymax, zmin, zmax)
--- GSL::Histogram3d#set_ranges_uniform([xmin, xmax], [ymin, ymax], [zmin, zmax])

--- GSL::Histogram3d#clone
--- GSL::Histogram3d#duplicate

--- GSL::Histogram3d#increment(x, y, z, weight = 1)
--- GSL::Histogram3d#accumulate(x, y, z, weight = 1)
--- GSL::Histogram3d#fill(x, y, z, weight = 1)
--- GSL::Histogram3d#increment2(x, y, z, weight = 1)
--- GSL::Histogram3d#accumulate2(x, y, z, weight = 1)
--- GSL::Histogram3d#fill2(x, y, z, weight = 1)
--- GSL::Histogram3d#get(i[, j[, k]])
--- GSL::Histogram3d#[]
--- GSL::Histogram3d#get_xrange(i)
--- GSL::Histogram3d#get_yrange(j)
--- GSL::Histogram3d#get_zrange(k)
--- GSL::Histogram3d#xrange
--- GSL::Histogram3d#yrange
--- GSL::Histogram3d#zrange
--- GSL::Histogram3d#bin
--- GSL::Histogram3d#xmax
--- GSL::Histogram3d#xmin
--- GSL::Histogram3d#ymax
--- GSL::Histogram3d#ymin
--- GSL::Histogram3d#zmax
--- GSL::Histogram3d#zmin
--- GSL::Histogram3d#nx
--- GSL::Histogram3d#ny
--- GSL::Histogram3d#nz
--- GSL::Histogram3d#reset
--- GSL::Histogram3d#find(x, y, z)
--- GSL::Histogram3d#max_val
--- GSL::Histogram3d#max_bin
--- GSL::Histogram3d#min_val
--- GSL::Histogram3d#min_bin
--- GSL::Histogram3d#xmean
--- GSL::Histogram3d#xsigma
--- GSL::Histogram3d#ymean
--- GSL::Histogram3d#ysigma
--- GSL::Histogram3d#zmean
--- GSL::Histogram3d#zsigma

--- GSL::Histogram3d#sum
--- GSL::Histogram3d#integral
--- GSL::Histogram3d#add(h2)
--- GSL::Histogram3d#sub(h2)
--- GSL::Histogram3d#mul(h2)
--- GSL::Histogram3d#div(h2)
--- GSL::Histogram3d#scale(val)
--- GSL::Histogram3d#shift(val)
--- GSL::Histogram3d#+(h2)
--- GSL::Histogram3d#-(h2)
--- GSL::Histogram3d#*(h2)
--- GSL::Histogram3d#/(h2)
--- GSL::Histogram3d#fwrite(io)
--- GSL::Histogram3d#fwrite(filename)
--- GSL::Histogram3d#fread(io)
--- GSL::Histogram3d#fread(filename)

--- GSL::Histogram3d#xyproject(kstart = 0, kend = nz-1)
    Creates a (({GSL::Histogram2d})) object by projecting the 3D histogram 
    ((|self|)) onto the xy-plane over the z-range from ((|kstart|)) to ((|kend|)).

--- GSL::Histogram3d#xzproject(jstart = 0, jend = ny-1)
    Creates a (({GSL::Histogram2d})) object by projecting the 3D histogram 
    ((|self|)) onto the xz-plane over the y-range from ((|jstart|)) to ((|jend|)).

--- GSL::Histogram3d#yzproject(istart = 0, iend = nx-1)
    Creates a (({GSL::Histogram2d})) object by projecting the 3D histogram 
    ((|self|)) onto the yz-plane over the x-range from ((|istart|)) to ((|iend|)).

((<prev|URL:hist2d.html>))
((<next|URL:ntuple.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
