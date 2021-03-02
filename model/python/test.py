from msc_gcv import *

p = Param()

p.kb   = 1e-4/yrsec
p.tauc = 120
p.Cw   = 2
p.U    = 4.9/mmyr
p.a    = 2.0
p.L    = 100e3

plot(p)

show()
