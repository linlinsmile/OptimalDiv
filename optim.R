########################
# Author: Linlin Tian
########################

# TODO: Add a summary description of what problem 
#       we are trying to solve here
# We are solving the optimal dividend problem of compound Poisson model. We approximate the continuous time model by a controlled 
# Markov chain. We derive the discrete time dynamic programming principle. From the discrete dynamic programming principle, we use value 
#iteration to get the result V^h. If h convergence to 0, then V^h convergence to the result of original model.
#We are trying to use the value iteration. 

########################
# set up model parameters
########################

# TODO: Add a human readable description
#       for each parameter
nSteps = 1000
g      = 0.0000000001 #tolerence
c      = 0.01 #discount factor
B      = 5000 #the upper bound of x
h      = 0.1   # step size

p      = 5 #premuim rate
lambda = 3 #intensity rate
alpha  = 2 #claim size intensity

echp = exp(-c * h / p)
lhp  = lambda * h / p
lhp1 = 1 - lhp
eah  = exp(-alpha * h)
lph2=echp*lhp1
lph3=lhp*echp*(1-exp(-alpha*h))
x    = seq(from = 0, to = B, by = h) #discret [0,B] by length of h

# Why 300 here? Give it a meaningful
# name and set it at the top,
# with the other defined constants
Vold = 300 * rep(1, length(x)) #300 is a random value,because the algorithm says for any given initial value,V_n will
#convergence to real value 

Vnew    = vector()

# Now we start the value iteration. I use length(Vold)-1 is because when we are calculating Vnew[k], we need the value of Vold[k+1]
for (t in 1:nSteps){ 
   Vnew[1] = lph2 * Vold[2]+lph3*Vold[1]
  
   # Explain what this inner loop is doing, from the second one, Vnew(x)=max{A,B}, A=(1-lambda*h/p)*e^(-ch/p)*Vold(x+h)+lhp*echp*Vsum,
         #B=Vold(x-h)+h
   for (k in 2:(length( Vold ) - 1)){
      m    = 0:(k-1)   
      Vsum = sum(Vold[k-m] * (1 - eah) *  eah^m)
      
      # What is this test about?
     # Because Vnew=max{A,B} if A>B, then choose A, else choose B.
      if((lhp1 * Vold[k + 1] + Vsum * lhp ) * echp > Vold[k - 1] + h){
         Vnew[k] = ( lhp1 * Vold[k + 1] + Vsum * lhp ) * echp
      }
      else{
         Vnew[k] = Vold[k - 1] + h
      }
   }

   Vold <- Vold[-length(Vold)]
   g[t] = max(abs(Vold - Vnew)[1:1000])
   Vold = Vnew
}
#here g[t] is to see the difference between Vold and Vnew, it is supposed to be very small, maybe less than 10^-6. In here I only 
#restrict the difference of Vold and Vnew for  the first 1000 number.
