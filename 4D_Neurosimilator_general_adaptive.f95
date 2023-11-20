program Fourth_order
! This program solves system of four ODE using 4th order Runge-Kutta method with variable step size
!  dv/dt = INPUT + (A - B)*v - (A+1)*w - B*y - H1 - H4
!  dw/dt = (A+1)*(v - w) - C*(exp(D*w) - 1.0)
!  dk/dt = H2 + H3 
!  dy/dt = G*(v - y)
!  To convert to 3rd order model, set H4=0 and dy/dt=0

implicit none
  
    ! Declare variables and parameters
    real, parameter :: Vth = 0.7, V_ref1 = 3.6, Vs = 10.0  ! fixed
    real, parameter :: t_end = 20.0                        
    real :: v, w, n, k, y, dt, t  
    real :: dv1, dv2, dv3, dv4, dw1, dw2, dw3, dw4
    real :: dk1, dk2, dk3, dk4, dy1, dy2, dy3, dy4           
    real :: H1, H2, H3, H4, Hs
    real :: INPUT    
    real, parameter :: A = 50                ! fixed
    real, parameter :: B = 0.16              ! 0.16=bursting/spiking, 0.12=chaotic 
    real, parameter :: C = 0.00005           ! fixed      
    real, parameter :: D = 20                ! fixed
    real, parameter :: E = 500               ! fixed
    real, parameter :: F = 32.08             ! 31.5=spiking; 32.08=bursting
    real, parameter :: G = 2.5               ! fixed    
    real, parameter :: V_ref2 = 0.85         ! fixed
    integer :: nsteps
    real :: max_deriv
    real :: precision, new_dt

    ! Set initial conditions
    v = 0.0
    w = 0.0
    k = 0.0
    y = 0.0
    dt = 0.01
    INPUT = 0.0
    
  
    ! Open output file for writing history data
    open(unit=1, file="history.dat", status="replace")
  
    ! Initialize time step and step counter
    t = 0.0
    nsteps = 0
    

    ! Time integration loop with variable step size using RK4 method
   
    do while (t < t_end)
      ! Calculate n value
      n = v - w

      ! Write current values 
      write (1,*) t, v, INPUT, w, n, k, y
  
      ! Set INPUT: for transient INPUT values it should be programmed apropriately
      if (t > 0.3) then
        INPUT = 0.28     
      else
        INPUT = 0
      endif  
      
      
      ! Calculate Heaviside functions H1, H2, H3 and H4
      if (k > Vth) then
        H1 = E*(k - Vth)
      else
        H1 = 0.0
      endif

      if (v >= V_ref1) then
        Hs = Vs
      else
        Hs = 0
      endif
      H2 = F*(Hs - k)
      
      
      if (k > Vth) then
        H3 = F*(Vth - k)
      else
        H3 = 0.0
      endif

      if (y >= V_ref2)  then     
        H4 = min(4.0, (E*(Vs - Vth)))      
      else 
        H4 = 0
      endif

      !if (v < 0 .and. ((dv1+dv2+dv3+dv4)/4)>0)  then
      !  H4 = -1.0*v                                   ! -1.9=bursting; -1.0=spiking;
      !else 
      !  H4 = H4
      !endif
      
        
      ! Calculate k1 values
      dv1 = INPUT + (A - B)*v - (A+1)*w - B*y - H1 - H4
      dw1 = (A+1)*(v - w) - C*(exp(D*w) - 1.0)
      dk1 = H2 + H3 
      dy1 = G*(v - y) 
      
      ! Calculate k2 values
      dv2 = INPUT + (A - B)*(v + 0.5*new_dt*dv1) - (A + 1)*(w + 0.5*new_dt*dw1) - B*(y + 0.5*new_dt*dy1) - H1 - H4
      dw2 = (A + 1)*((v + 0.5*new_dt*dv1) - (w + 0.5*new_dt*dw1)) - C*(exp(D*(w + 0.5*new_dt*dw1)) - 1.0)
      dk2 = H2 + H3 
      dy2 = G*((v + 0.5*new_dt*dv1) - (y + 0.5*new_dt*dy1)) 
      
      ! Calculate k3 values
      dv3 = INPUT + (A - B)*(v + 0.5*new_dt*dv2) - (A + 1)*(w + 0.5*new_dt*dw2) - B*(y + 0.5*new_dt*dy2) - H1 - H4
      dw3 = (A + 1)*((v + 0.5*new_dt*dv2) - (w + 0.5*new_dt*dw2)) - C*(exp(D*(w + 0.5*new_dt*dw2)) - 1.0)
      dk3 = H2 + H3 
      dy3 = G*((v + 0.5*new_dt*dv2) - (y + 0.5*new_dt*dy2)) 
      
      ! Calculate k4 values
      dv4 = INPUT + (A - B)*(v + new_dt*dv3) - (A + 1)*(w + new_dt*dw3) - B*(y + new_dt*dy3) - H1 - H4
      dw4 = (A + 1)*((v + new_dt*dv3) - (w + new_dt*dw3)) - C*(exp(D*(w + new_dt*dw3)) - 1.0)
      dk4 = H2 + H3 
      dy4 = G*((v + new_dt*dv3) - (y + new_dt*dy3)) 

      ! Calculate the maximum derivative among all variables:
      max_deriv = max(abs((dv1+dv2+dv3+dv4)/4), abs((dw1+dw2+dw3+dw4)/4), abs((dk1+dk2+dk3++dk4)/4), abs((dy1+dy2+dy3+dy4)/4))

      ! Calculate the variable time step size:
      !  step size is recalculated within each iteration of the time integration loop 
      !  based on the maximum derivative among all variables
      precision = 0.0001                      ! in Volts
      if (max_deriv > 0.0) then
         new_dt = precision / max_deriv       ! V/(dV/dt)= dt
      else
      new_dt = dt  ! Use a default time step size or skip step size update
      endif


     ! Update v, w, and k variables using k1/k2/k3/k4 values and Runge-Kutta formula
     v = v + (1.0/6.0)*(dv1 + 2.0*dv2 + 2.0*dv3 + dv4)*new_dt
     w = w + (1.0/6.0)*(dw1 + 2.0*dw2 + 2.0*dw3 + dw4)*new_dt
     k = k + (1.0/6.0)*(dk1 + 2.0*dk2 + 2.0*dk3 + dk4)*new_dt
     y = y + (1.0/6.0)*(dy1 + 2.0*dy2 + 2.0*dy3 + dy4)*new_dt

     ! Update time and step counter
     t = t + new_dt
     nsteps = nsteps + 1

    end do

! Close output file
close(1)

! Output final values of variables
print *, "Final values of v, n, k, y, INPUT:"
print *, v, n, k, y, INPUT
print *, "Number of time steps:", nsteps

end program Fourth_order