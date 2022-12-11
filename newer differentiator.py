"""
1d-time independent schrodinger equation plotter
Sean Heffley 2022
PHY446
"""

# a library for finding eigenstates/valid energy levels : https://github.com/FelixDesrochers/Numerov/blob/master/Fct_Numerov.py
# has a lot of functions verifying conditions that would be helpful expending this to handle an arbitrary potential function and automatically plot it but assumes mass = 1 (electron)
#IMPORTS/////////////////////////////////////////////////////////////////////////////////////
#I like this because it leaves a readable message if run outside an ide (command prompt) and doesn't require popups
try:import numpy as np
except:input("failed to import numpy")
try:
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button
except:input("failed to import matplotlib")
try:from scipy.integrate import odeint,simps
except:input("failed to import scipy")
try:from decimal import Decimal
except:input("failed to import decimal")

#initiate INDEPENDENTS////////////////////////////////////////////////////////////////////////////////
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

#ivp
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

#PLOTTING/////////////////////////////////////////////////////////////////////////////////////////////// warning messy spagetti code ahead ////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    
    #setup##########################################################################################
    fig, ax = plt.subplots(figsize = (12,8))
    arrayX = np.arange(x_min,x_max,.02)
    plt.xlabel('x')
    
    #plotting graphs##############################################################
    arrayOdePsi = waveFunctionSystemODE(W0, arrayX)#and global variables mass and energy   
    arrayNumPsi = waveFunctionNumerov(augment, arrayX)#and global variables mass and energy
    
    potential_plot    = ax.plot(arrayX, arrayFunction(arrayX,potential), lw=2, label='Potential', color = 'black')    
    wave_odeint_plot  = ax.plot(arrayX, arrayOdePsi,  label='ODEint solution wave function', color='cyan')                      
    wave_numerov_plot = ax.plot(arrayX, arrayNumPsi,  label='Numerov solution wave function',color='orange')
    
    ax.set_ylim(-4,6)
    leg = ax.legend()

    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.5)

    #integrate over solution^2 and display sqrt(of the area) --- normalization for wave funciton
    odeintegral = simps([y**2 for y in arrayOdePsi],arrayX)
    axodeintegral = fig.add_axes([0.69,0.12,0.1,0.03])
    odeintegralbutton = Button(axodeintegral, '%.2E' % Decimal(str(np.sqrt(odeintegral))), color = "cyan",hovercolor = "cyan")

    numintegral = simps([y**2 for y in arrayNumPsi],arrayX)
    axnumintegral = fig.add_axes([0.58,0.12,0.1,0.03])
    numintegralbutton = Button(axnumintegral, '%.2E' % Decimal(str(np.sqrt(numintegral))), color = "orange",hovercolor = "orange")


    #dynamic graph stuff#################################################### this is where the mess begins (first time using matplotlib beyond x,y arrays autoformat)
    #making space
    fig.subplots_adjust(left=.14,right = 0.8, bottom=0.2, top = .9)

    #toggle solutions with normalizations---------------------------
    global plot_ode
    global plot_num #this section would be abomidable as a library wouldn't it
    plot_ode = True
    plot_num = True

    def toggle_ode(event):
        global plot_ode
        plot_ode = not plot_ode
        update(None)
        
    def toggle_num(event):
        global plot_num
        plot_num = not plot_num
        update(None)
        
    odeintegralbutton.on_clicked(toggle_ode)
    numintegralbutton.on_clicked(toggle_num)
    
    #toggle normalize button---------------------------- (lambda and a make toggle button function perhaps?)
    global normalize
    normalize = False #again, usure how to do this without global variables

    def toggle_normalize(event):
        global normalize
                
        if normalizebutton.color == "green":
            normalizebutton.color = "red"
            normalize = False
        else:
            normalizebutton.color = "green"
            normalize = True

        update(None)
            
    axnormalize = fig.add_axes([0.15,0.12,0.1,0.03])
    normalizebutton = Button(axnormalize, "toggle normalize", hovercolor = '0.975', color = "red")
    normalizebutton.on_clicked(toggle_normalize)

    #toggle limit axes button----------------------------
    global limit_yaxis
    limit_yaxis = True 

    def toggle_axis_limit(event):
        global limit_yaxis
                
        if axisbutton.color == "green":
            axisbutton.color = "red"
            limit_yaxis = False
        else:
            axisbutton.color = "green"
            limit_yaxis = True
            
        update(None)
    
    axlimit = fig.add_axes([0.26,0.12,0.1,0.03])
    axisbutton = Button(axlimit, "toggle axis_limit", hovercolor = '0.975', color = "green")
    axisbutton.on_clicked(toggle_axis_limit)

    #plot probaility toggle---------------------------
    global plot_probability
    plot_probability = False
    
    def toggle_plot_probability(event):
        global plot_probability
                
        if plot_probabilitybutton.color == "green":
            plot_probabilitybutton.color = "red"
            plot_probability = False
        else:
            plot_probabilitybutton.color = "green"
            plot_probability = True

        update(None)
            
    axplot_probability = fig.add_axes([0.81,0.12,0.15,0.03])
    plot_probabilitybutton = Button(axplot_probability, "toggle plot_probability", hovercolor = '0.975', color = "red")
    plot_probabilitybutton.on_clicked(toggle_plot_probability)


    #sliders-----------------------------------------    
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
        color = 'cyan',
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
        color = 'cyan',
    )

    def update(val):######################################################################
        #clear plot
        ax.clear()
        ax.set_xlabel('x')

        #update values
        global energy
        global mass
        energy = energy_slider.val
        mass = mass_slider.val

        potential_plot    = ax.plot(arrayX, arrayFunction(arrayX,potential), lw=2, label='Potential', color = 'black')
        
        if plot_num:######################################################################  
            arrayNumPsi = waveFunctionNumerov(augment_slider.val, arrayX)#and global variables mass and energy

            #update normalization---------------------------
            numintegral = simps([y**2 for y in arrayNumPsi],arrayX)
            numintegralbutton.label.set_text('%.2E' % Decimal(str(np.sqrt(numintegral))))

            if normalize:
                #divide each element of wave solutions by normalization
                arrayNumPsi[:] = [y/np.sqrt(numintegral) for y in arrayNumPsi]

                #sanity check
                numintegral = simps([y**2 for y in arrayNumPsi],arrayX)
                numintegralbutton.label.set_text('%.2E' % Decimal(str(np.sqrt(numintegral))))

            #re-plot wave function                      
            wave_numerov_plot = ax.plot(arrayX, arrayNumPsi,  label='Numerov solution wave function',color='orange')

            if plot_probability:                  
                wavesquared_numerov_plot = ax.plot(arrayX, [y**2 for y in arrayNumPsi],  label='Numerov solution probaility density',color='red')

        if plot_ode:######################################################################
            arrayOdePsi = waveFunctionSystemODE([psi0_slider.val,dpsi0_slider.val], arrayX)#and global variables mass and energy   

            #update normalization----------------------------
            odeintegral = simps([y**2 for y in arrayOdePsi],arrayX)
            odeintegralbutton.label.set_text('%.2E' % Decimal(str(np.sqrt(odeintegral))))

            if normalize:
                #divide each element of wave solutions by normalization
                arrayOdePsi[:] = [y/np.sqrt(odeintegral) for y in arrayOdePsi]

                #sanity check
                odeintegral = simps([y**2 for y in arrayOdePsi],arrayX)
                odeintegralbutton.label.set_text('%.2E' % Decimal(str(np.sqrt(odeintegral))))

            #re-plot wave function    
            wave_odeint_plot  = ax.plot(arrayX, arrayOdePsi,  label='ODEint solution wave function', color='cyan')                      

            if plot_probability:
                wavesquared_odeint_plot  = ax.plot(arrayX, [y**2 for y in arrayOdePsi],  label='ODEint solution probability density', color='blue')                      

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


    #searching energy values manually. library noted at top does automatically but limits mass = 1
    #set min/max/reset (no undo)

##    energy_slider = Slider(
##        ax=axenergy,
##        label='energy',
##        valmin=0,
##        valmax=2,
##        valinit=energy,
##        color = 'yellow',

    def set_energy_max(event):#namespace is getting disgusting lol, also not optimized at all. would have to do a lot more work with matplotlib to comfortably do that here
        energy_slider.valmax = energy_slider.val
        energy_slider.ax.set_xlim(energy_slider.valmin,energy_slider.valmax)
        update(None)
        
    def set_energy_min(event):
        energy_slider.valmin = energy_slider.val
        energy_slider.ax.set_xlim(energy_slider.valmin,energy_slider.valmax)
        update(None)

    def reset_energy_bounds(event):
        energy_slider.valmin = 0
        energy_slider.valmax = 2
        energy_slider.ax.set_xlim(energy_slider.valmin,energy_slider.valmax)
        update(None)
        
    axset_energy_max = fig.add_axes([0.1,0.02,0.12,0.03])
    axset_energy_min = fig.add_axes([0.25,0.02,0.12,0.03])
    axreset_energy_bounds = fig.add_axes([0.4,0.02,0.12,0.03])

    set_energy_maxbutton = Button(axset_energy_max, "set_energy_max", hovercolor = '0.975', color = "grey")
    set_energy_minbutton = Button(axset_energy_min, "set_energy_min", hovercolor = '0.975', color = "grey")
    reset_energy_boundsbutton = Button(axreset_energy_bounds, "reset energy slider", hovercolor = '0.975', color = "grey")

    set_energy_maxbutton.on_clicked(set_energy_max)
    set_energy_minbutton.on_clicked(set_energy_min)
    reset_energy_boundsbutton.on_clicked(reset_energy_bounds)
    
    plt.show()
    
