#!/usr/bin/env ruby
require("gsl")

n = 256
nc = 20

data = GSL::Vector.alloc(n)
data.fscanf("ecg.dat")

w = GSL::Wavelet.alloc("daubechies", 4)
work = GSL::Wavelet::Workspace.alloc(n)

# Choose as you like...
data2 = w.transform(data, GSL::Wavelet::FORWARD, work)
#data2 = data.wavelet_transform(w, GSL::Wavelet::FORWARD, work)
#data2 = data.wavelet_transform_forward(w, work)
#data2 = w.transform(data, work)
#data2 = w.transform(data)
#data2 = w.transform(data, GSL::Wavelet::FORWARD)
#data2 = w.transform_forward(data, work)
#data2 = w.transform_forward(data)
#data2 = GSL::Wavelet.transform(w, data, GSL::Wavelet::FORWARD, work)
#data2 = GSL::Wavelet.transform(w, data, GSL::Wavelet::FORWARD)
#data2 = GSL::Wavelet.transform(w, data, work)
#data2 = GSL::Wavelet.transform(w, data)
#data2 = GSL::Wavelet.transform_forward(w, data, work)
#data2 = GSL::Wavelet.transform_forward(w, data)

perm = data2.abs.sort_index

i = 0
while (i + nc) < n
  data2[perm[i]] = 0.0
  i += 1
end

# Choose as you like...
data3 = w.transform(data2, GSL::Wavelet::BACKWARD, work)
#data3 = data2.wavelet_transform(w, GSL::Wavelet::BACKWARD, work)
#data3 = data2.wavelet_transform_inverse(w, work)
#data3 = w.transform(data2, GSL::Wavelet::BACKWARD)
#data3 = w.transform_inverse(data2, work)
#data3 = w.transform_inverse(data2)
#data3 = GSL::Wavelet.transform(w, data2, GSL::Wavelet::BACKWARD, work)
#data3 = GSL::Wavelet.transform(w, data2, GSL::Wavelet::BACKWARD)
#data3 = GSL::Wavelet.transform_inverse(w, data2, work)
#data3 = GSL::Wavelet.transform_inverse(w, data2)

GSL::graph(nil, data, data3, "-T X -C -g 3 -x 0 #{data.size} -L 'Red: data, Green: DWT'")

