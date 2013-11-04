      function bi1(z)
      include 'mach.p'
      real*8 bi1,z,t
c modified bessel function I_1(z)/z
      t=.071111111111111*z**2
      bi1=.5+t*(.87890594+t*(.51498869+t*(.15084934+t*(.02658733
     ++t*(.00301532+t*.00032411)))))
      return
      end
