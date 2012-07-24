#!/usr/bin/env ruby
require("gsl")

data =  GSL::Vector[17.2, 18.1, 16.5, 18.3, 12.6]
mean     = data.mean()
variance = data.variance()
largest  = data.max()
smallest = data.min()

printf("The dataset is %g, %g, %g, %g, %g\n",
        data[0], data[1], data[2], data[3], data[4]);

printf("The sample mean is %g\n", mean);
printf("The estimated variance is %g\n", variance);
printf("The largest value is %g\n", largest);
printf("The smallest value is %g\n", smallest);

p data.stats_sd
