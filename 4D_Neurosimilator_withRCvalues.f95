program fourD_Neurosimilator

! This program solves system of four ODE using 4th order Runge-Kutta method
!  v' = INPUT + (A-B)*v - (A+1)*w - B*y - H1 - H4
!  w' = (A+1)*(v - w) - C*(exp(D*w) - 1.0)
!  k' = H2 + H3 
!  y' = G*(v - y)
!  where:
!  H1 = E*(k-Vth)*heav(k-Vth)
!  H2 = F*(heav(v-Vref1)*Vs-k)
!  H3 = F*(Vth-k)*heav(k-Vth)
!  H4 = min(4.0, (E*(Vs-Vth))*heav(y-Vref2)) 
!    
!  With the constants:
!  INPUT = Vin/(Rin*Cm)
!  A = 1/(Rg*Cm) - 1/(Rm*Cm) - 1/(Rin*Cm)   =>+/- 50.0
!  A+1 = 1/(Rg*Cm) = +/- = 1/(Rg*C1)        =>+/- 51.0
!  B = 1/(R3*Cm)                            =>+/- 0.16 for spiking/bursting and 0.12 for chaotic 
!  C = Is/C1                                =>0.00005
!  D = 1.6E-19/(1.9*1.38E-23*300)           =>20.0
!  E = hfe/(R2*Cm)= hfe/(R4*Cm)             =>500
!  F = 1/(R1*C2)                            =>31.5 for spiking and 32.08 for bursting/chaotic
!  G = 1/(R3*C3)                            =>2.5

  implicit none
  
    ! Declare variables and parameters
    real, parameter :: Vth = 0.7, V_ref1 = 3.6, V_ref2 = 0.9, Vs = 10.0  ! threshold/reference/positive rail voltages
    real, parameter :: dt = 5.0E-6, t_end = 5.0                          ! integration parameters
    real :: v               ! membrane potential dimension variable
    real :: w, n            ! depolarization dimensions: n=(v-w) mimics m^3*h in HH model
    real :: k               ! fast repolarization dimension: k mimics n^4 in HH model
    real :: y               ! slow repolarization dimension
    real :: t, H1, H2, H3, H4, Hs!, Vin                                 ! time, Heavisides and input voltage variables
    real :: dv1, dv2, dv3, dv4, dw1, dw2, dw3, dw4                      ! rk4 values
    real :: dk1, dk2, dk3, dk4, dy1, dy2, dy3, dy4                      ! rk4 values
    real :: INPUT, A, B, C, D, E, F, G                                  ! diff equations constants
    real, parameter :: Cm = 19.6E-6                                     ! membrane capacitance (v-function)
    real, parameter :: C1 = 20.0E-6                                     ! w/n-function capacitance
    real, parameter :: C2 = 3.3E-6                                      ! k-function capacitance
    real, parameter :: R1 = 9.63E3, R2 = 10.0E3                         ! k-function resistances (in case of spiking R1=9.45E3)
    real, parameter :: Rin = 100.0E3                                    ! input resistance
    real, parameter :: Rm = 100.0E3                                     ! membrane resistance
    real, parameter :: Rg = 1000                                        ! w/n-function resistance
    real, parameter :: R3 = 320.0E3, R4 = 10.0E3                        ! y function resistances
    real, parameter :: C3 = 1.25E-6                                     ! y function capacitance
    real, parameter :: Is = 1.0E-9                                      ! diode saturation current
    real, parameter :: hfe = 100                                        ! transistor current gain
    integer :: nsteps
    real :: new_dt !, max_deriv, precision                              ! activate "max_deriv" and "precision" for adaptive step calculations

    ! Set initial conditions
    v = 0.0
    w = 0.0
    k = 0.0
    y = 0.0
 !  dt = 0.01         ! activate for adaptive time step                 
    new_dt = dt       ! inactivate for adaptive  time step
    INPUT = 0.0

    ! Open output file for writing history data
    open(unit=1, file="history.dat", status="replace")
  
    ! Initialize time step and step counter
    t = 0.0
    nsteps = 0

    ! Calculate constants
    A = 1/(Rg*Cm) - 1/(Rm*Cm) - 1/(Rin*Cm)  
    ! A+1 = 1/(Rg*Cm) = +/- = 1/(Rg*C1)        
    B = 1/(R3*Cm)                             
    C = Is/C1                               
    D = 1.6E-19/(1.9*1.38E-23*300)           
    E = hfe/(R2*Cm) ! = hfe/(R4*Cm)             
    F = 1/(R1*C2)                            
    G = 1/(R3*C3)
  
       ! Time integration loop with variable step size using RK4 method
   
    do while (t < t_end)
        ! Calculate n value
        n = v - w
  
        ! Write current values 
        write (1,*) t, v, INPUT, w, n, k, y
    
        ! Set INPUT: for transient INPUT values it should be programmed apropriately
        if (t > 0.3) then
          INPUT = 0.5     
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

    ! activate the code below for adaptive time step calculations
  
        ! Calculate the maximum derivative among all variables:
    !    max_deriv = max(abs((dv1+dv2+dv3+dv4)/4), abs((dw1+dw2+dw3+dw4)/4), abs((dk1+dk2+dk3++dk4)/4), abs((dy1+dy2+dy3+dy4)/4))
  
        ! Calculate the variable time step size:
        !  step size is recalculated within each iteration of the time integration loop 
        !  based on the maximum derivative among all variables

    !    precision = 0.0001                      ! in Volts
    !    if (max_deriv > 0.0) then
    !       new_dt = precision / max_deriv       ! V/(dV/dt)= dt
    !    else
    !    new_dt = dt  ! Use a default time step size or skip step size update
    !    endif
  
  
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
    print *, "Final values of v, n, k, y, INPUT, A, B, C, D, E, F, G:"
    print *, v, n, k, y, INPUT, A, B, C, D, E, F, G
    print *, "Number of time steps:", nsteps
    
    end program fourD_Neurosimilator