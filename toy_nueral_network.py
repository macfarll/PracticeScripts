'''
Python implementation of a toy neural network, as inspired by:
http://karpathy.github.io/neuralnets/
'''

#Import math and create a sigmoid function for future use
import math
def sig(x):
    return 1/(1+math.exp(-x))
#step size will be used later to determine the rate at which
# the network modifies the input
step_size=0.01

#Define Units, useful for passing forward VALues and backprop GRADients
class Unit:
    def __init__(self,val,grad):
        self.val=val
        self.grad=grad

#A template for an addition gate
class AddGate:
#The forward calculation where the sum is saved to top.val
    def forward(self,x,y):
        self.x=x
        self.y=y
        top=x.val+y.val
        self.top=Unit(top,0)
        return self.top
    def backward(self):
#The backprop where the gradient in respect to x is saved in x.grad
# and the gradient with respect to y in in y.grad, these are dependent
# on only the input gradient, since the inputs for the gate affect the
# output in a linear fashion, the derivative of each input with 
# respect to the output is 1.
        self.x.grad+=self.top.grad
        self.y.grad+=self.top.grad

#A template for forward multiplication gates
class MultiGate:
#The forward calculation, where the product on the inputs is in top.val
    def forward(self,x,y):
        self.x=x
        self.y=y
        top=self.x.val*self.y.val
        self.top=Unit(top,0)
        return self.top
#Backpropogation where the change in the gate output relative to
# each input is stored in input.grad. The derivative of the input is
# represented by the (output with respect to the final product, 
# top.grad) multiplied by (the input in respect to the output, which
# for x would be y, and for y would be x) 
    def backward(self):
        self.x.grad+= self.y.val*self.top.grad
        self.y.grad+= self.x.val*self.top.grad

#A gate for a sigmoid function
class SigGate():
    def forward(self,x):
        self.x=x
        top=sig(x.val)
        self.top=Unit(top,0)
        return self.top
#A bit more complex calculus, but still just the output with respect
# to the input
    def backward(self):
        s=sig(self.x.val)
        self.x.grad+=(s*(1-s))*self.top.grad


#Starting values
u=list()
u.append(Unit(1,0))
u.append(Unit(3,0))
u.append(Unit(-4,0))
u.append(Unit(-1,0))

#Starting gates
mulg1=MultiGate()
mulg2=MultiGate()
addg1=AddGate()
sigg=SigGate()

#Function to run the circuit forward and return the final output
# Function is: ax + by + c
# or....       
def fwdNeuron():
    passInstance=list()
    XnY= mulg1.forward(u[0],u[1])
    XnYpA= addg1.forward(XnY,u[2])
    PnZ= mulg2.forward(XnYpA,u[3])
    SP= sigg.forward(PnZ)
    print 'Output for current pass: ', SP.val
    return SP

#Function to start the backprop
def backNeuron():
    sigg.backward()
    mulg2.backward()
    addg1.backward()
    mulg1.backward()

#Finally, we run the network. First, we proceed forward through,
# next we initialize the backprop*, and finally we iterate through
# our list of inputs and adjust them slightly. The output will only
# be the output of the forward gradient, a subsequent run of main()
# will return the adjusted output. 

# *-the backprop is initialized at 1 since the final ouput of the 
# network with respect to the output of the final gate is 1, since
# they are the same number.  
def __main__():
    endNode=fwdNeuron()
    endNode.grad=1
    backNeuron()
    for i in u:
        i.val +=step_size*i.grad

#Run main 10 times
for i in range(10):
    __main__()
