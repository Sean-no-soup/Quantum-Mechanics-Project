"""
1d-time independent schrodinger equation plotter
Sean Heffley 2022
PHY446
"""

# a library for finding eigenstates/valid energy levels : https://github.com/FelixDesrochers/Numerov/blob/master/Fct_Numerov.py
# has a lot of functions verifying conditions that would be helpful expending this to handle an arbitrary potential function and automatically plot it
#IMPORTS/////////////////////////////////////////////////////////////////////////////////////
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
#import matplotlib.cm as cm
from scipy.integrate import odeint

#INDEPENDENTS////////////////////////////////////////////////////////////////////////////////
#global variables probably aren't the most elegant solution but having trouble tracking down how to avoid them in the scipy implementation
global energy
global mass

energy = 0.1
mass = 36

#doesn't check for conditions just graphs based on these,
#i.e. automatically determining a reasonable window for the function and initial values not implemented
x_min = -4 
x_max = 4


def potential(x): #u0  *  (1-np.e **  (-(x**2)/(a**2)))
    u0 = 1
    a = 1
    return (u0*(1 - np.e ** (  (-(x)**2)  /  ((a)**2)  )))

#initial values
W0 = [0,0.5] #scipy diffeq
augment = 0.01 #numerov 0.01 vs [0,.5] are close 


#GENERAL FUNCTIONS////////////////////////////////////////////////////////////////////////////
def arrayFunction(arrayX, F):
    arrayF = []
    for i in arrayX:
        arrayF.append(F(i)) 
    return arrayF

def g(u): #2m(e-u) ,aka k**2, used in both methods
    return 2*mass*(energy-u)


#NUMEROV METHOD///////////////////////////////////////////////////////////////////////////////
def waveFunctionNumerov(augment, arrayX):
    """returns a list of psi (floats)\

    parameters\
        augment : delta psi between the first 2 points. if exactly == 0 psi remains 0 over any domain?!!!
        arrayX : evenly spaced array of x values associated with psi

        energy global var
        mass global var
        
    potential(x) function defined elsewhere 
    """
    waveFunction = [] #making a list for psi values

    #need spacing
    xs = arrayX[1]-arrayX[0]
    
    #setting x(i), x(i+1)
    waveFunction.append(0)
    waveFunction.append(augment)

    #LOOOOOP
    for index in range(len(arrayX)-1): #truncates last point to avoid index errors
        if index > 0: #skip first , second adds third point, ect. 
            
            X_plus  = arrayX[index+1]
            X      = arrayX[index]
            X_minus = arrayX[index-1]
            
            #potential
            U_plus , U, U_minus = map(potential, [X_plus,X,X_minus])

            #2m(U-E), aka k**2
            G_plus , G, G_minus = map(g, [U_plus,U,U])

            #psi_plus                              psi                                            psi_minus
            psi_plus = ((2*(1-(5/12)*(xs**2)*(G))*(waveFunction[-1]))-(1+(1/12)*(xs**2)*G_minus)*(waveFunction[-2]))/(1+(1/12)*(xs**2)*G_plus)

            #add next point
            waveFunction.append(psi_plus)

    return waveFunction



#FIRST ORDER SYSTEM METHOD///////////////////////////////////////////////////////////////////////////
#[dw0dx] == [w1]
#[dw1dx] == [w2] == -(2m/(hbar^2))*(e-u)*[w0] == -g[w0]

#vector system
def dWdx(W, x):
    dw0dx = W[1]
    dw1dx = -g(potential(x)) * W[0]
    return [dw0dx, dw1dx]

def waveFunctionSystemODE(W0, arrayX):
    """returns a python list of psi (floats)\
    parameters\
    
        arrayX: iterable x values associated with psi
        W0 : initial wave value and derivative in a list [a,b]

        global mass var
        global energy var
        
    potential(x) function defined elsewhere 
    """
    #solve system #from scipy.org# scipy.integrate.odeint(func, y0, t, args=(), Dfun=None, col_deriv=0, full_output=0, ml=None, mu=None, rtol=None, atol=None, tcrit=None, h0=0.0, hmax=0.0, hmin=0.0, ixpr=0, mxstep=0, mxhnil=0, mxordn=12, mxords=5, printmessg=0, tfirst=False)
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
    W = odeint(dWdx, W0, arrayX)    #returns an array, values  psi and dpsi/dx per entry
    
    return list(W[:,0]) #python list of the first column

#PLOTTING///////////////////////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    
    #setup##########################################################################################
    fig, ax = plt.subplots(figsize = (12,8))
    arrayX = np.arange(x_min,x_max,.02)
    plt.xlabel('x')

    
    #plotting graphs##############################################################
    potential_plot    = ax.plot(arrayX, arrayFunction(arrayX,potential), lw=2, label='Potential', color = 'black')
    wave_odeint_plot  = ax.plot(arrayX, waveFunctionSystemODE(W0, arrayX),     label='ODEint solution', color='blue') #and global variables mass and energy                        
    wave_numerov_plot = ax.plot(arrayX, waveFunctionNumerov(augment, arrayX),  label='Numerov solution',color='orange') #and global variables mass and energy

    ax.set_ylim(-4,6)
    leg = ax.legend()

    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.5)

    #dynamic graph stuff####################################################
    #toggle limit axes button. for dynamic graph####################################################
    global limit_yaxis
    limit_yaxis = True #again, usure how to do this without global variables

    def toggle_axis_limit(event):
        global limit_yaxis
                
        if button.color == "green":
            button.color = "red"
            limit_yaxis = False
        else:
            button.color = "green"
            limit_yaxis = True
    
    axlimit = fig.add_axes([0.7,0.12,0.1,0.03])
    button = Button(axlimit, "toggle axis_limit", hovercolor = '0.975', color = "green")
    button.on_clicked(toggle_axis_limit)

    #sliders#############################################################    
    #making space
    fig.subplots_adjust(left=.14,right = 0.8, bottom=0.2, top = .9)
    
    #energy
    axenergy = fig.add_axes([.1, 0.07, 0.8, 0.03])
    energy_slider = Slider(
        ax=axenergy,
        label='energy',
        valmin=0,
        valmax=2,
        valinit=energy,
        color = 'yellow',
    )

    #mass
    axmass = fig.add_axes([0.04, 0.2, 0.03, 0.7])
    mass_slider = Slider(
        ax=axmass,
        label="mass",
        valmin=0,
        valmax=40,
        valinit=mass,
        orientation="vertical",
        color = 'grey',
    )

    #augment
    axaugment = fig.add_axes([0.82, 0.2, 0.03, 0.7])
    augment_slider = Slider(
        ax=axaugment,
        label="psi(i=1)",
        valmin=-1,
        valmax=1,
        valinit=augment,
        orientation="vertical",
        color = 'orange',
    )

    #psi0
    axpsi0 = fig.add_axes([0.87, 0.2, 0.03, 0.7])
    psi0_slider = Slider(
        ax=axpsi0,
        label="psi(i=0)",
        valmin=-1,
        valmax=1,
        valinit=W0[0],
        orientation="vertical",
        color = 'blue',
    )

    #dpsi0
    axdpsi0 = fig.add_axes([0.92, 0.2, 0.03, 0.7])
    dpsi0_slider = Slider(
        ax=axdpsi0,
        label="dpsi(i=0)",
        valmin=-1,
        valmax=1,
        valinit=W0[1],
        orientation="vertical",
        color = 'blue',
    )

    def update(val):
        global energy
        global mass
        
        energy = energy_slider.val
        mass = mass_slider.val

        ax.clear()

        potential_plot    = ax.plot(arrayX, arrayFunction(arrayX,potential), lw=2, label='Potential', color = 'black')
        wave_odeint_plot  = ax.plot(arrayX, waveFunctionSystemODE([psi0_slider.val, dpsi0_slider.val], arrayX),     label='ODEint solution', color='blue') #and global variables mass and energy                        
        wave_numerov_plot = ax.plot(arrayX, waveFunctionNumerov(augment_slider.val, arrayX),  label='Numerov solution',color='orange') #and global variables mass and energy

        leg = ax.legend()

        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.5)

        if limit_yaxis: ax.set_ylim(-4,6)
        
        fig.canvas.draw_idle()

    # register the update function with each slider
    mass_slider.on_changed(update)
    energy_slider.on_changed(update)
    augment_slider.on_changed(update)
    psi0_slider.on_changed(update)
    dpsi0_slider.on_changed(update)

    """

    #manually search for energy levels (psi -> 0 as x-> +- infinity) energy, set min/max energy values of slider
    def set_energy_max(event):
    def set_energy_min(event):

    
    axlimit = fig.add_axes([0.7,0.12,0.1,0.03])
    button = Button(axlimit, "toggle axis_limit", hovercolor = '0.975', color = "green")
    button.on_clicked(toggle_axis_limit)
    """

    plt.show()
    
