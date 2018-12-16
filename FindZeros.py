import tkinter as tk
import tkinter.scrolledtext as tkst
import math
from cmath import log10, log, pi, e, sin, cos, tan, acos, asin, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp, sqrt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

userFunction=''

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.entry_text = tk.Label(self)
        self.entry_text["text"] = "f(x) = "
        self.entry_text.grid(row=0,column=0, padx=[20,0], pady=[20,0],sticky="W")

        self.entry_box = tk.Entry(self,width=19)
        self.entry_box.grid(row=0,column=1,columnspan=2,pady=[20,0])
        self.entry_box.insert(1,"x-3*x**4-x**6")

        self.tol_text = tk.Label(self)
        self.tol_text["text"] = "Tolerance "
        self.tol_text.grid(row=1,column=0,padx=[20,0],sticky="W")

        self.tol_box = tk.Entry(self,width=19)
        self.tol_box.grid(row=1,column=1,columnspan=2)
        self.tol_box.insert(1,"0.01")

        self.xlim_text = tk.Label(self)
        self.xlim_text["text"] = "Real Axis: "
        self.xlim_text.grid(row=2,column=0,padx=[20,0],sticky="W")
        
        self.xlim_M_box = tk.Entry(self,width=9)
        self.xlim_M_box.grid(row=2,column=1,padx=[20,0])
        self.xlim_M_box.insert(1,"-2")
        
        self.xlim_m_box = tk.Entry(self,width=9)
        self.xlim_m_box.grid(row=2,column=2,padx=[0,20])
        self.xlim_m_box.insert(1,"2")
        
        self.ylim_text = tk.Label(self)
        self.ylim_text["text"] = "Imaginary Axis: "
        self.ylim_text.grid(row=3,column=0,padx=[20,0],sticky="W")
        
        self.ylim_M_box = tk.Entry(self,width=9)
        self.ylim_M_box.grid(row=3,column=1,padx=[20,0])
        self.ylim_M_box.insert(1,"-2")
        
        self.ylim_m_box = tk.Entry(self,width=9)
        self.ylim_m_box.grid(row=3,column=2,padx=[0,20])
        self.ylim_m_box.insert(1,"2")

        self.error_text = tk.Label(self, fg="red")
        self.error_text["text"] = ""
        self.error_text.grid(row=4,column=0, columnspan=2, padx=[20,0])

        self.compute_button = tk.Button(self)
        self.compute_button.grid(row=4,column=2, padx=[0,20], pady=[5,5])
        self.compute_button["text"] = "Compute"
        self.compute_button["command"] = self.CollectAndCall
        
        self.output_area = tkst.ScrolledText(self,wrap=tk.WORD,width=15,height=7)
        self.output_area.grid(row=0,column=3, rowspan=5, padx=[0,20])
        self.output_area.insert(tk.INSERT,"(-0.38-0.59j)\n0j\n(-0.38+0.59j)\n(0.66+0j)\n(0.06+1.74j)\n(0.06-1.74j)")
        
        fig1 = Figure(figsize=(5,5), dpi=70)
        fig = fig1.add_subplot(111)
        fig.scatter([-0.38,0,-0.38,0.66,0.06,0.06],[0.59,0,-0.59,0,1.74,-1.74])

        self.canvas = FigureCanvasTkAgg(fig1,self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=6, columnspan=4, pady=[0,40])
        
    def CollectAndCall(self):
        try:
            tol=float(self.tol_box.get())
        except:
            self.error_text["text"] = "tol not a number"
            return
        try:
            a=float(self.xlim_M_box.get())
        except:
            self.error_text["text"] = "Real Max Bound not a number"
            return
        try:
            b=float(self.xlim_m_box.get())
        except:
            self.error_text["text"] = "Real Min Bound not a number"
            return
        try:
            c=float(self.ylim_M_box.get())
        except:
            self.error_text["text"] = "Imaginary Max Bound not a number"
            return
        try:
            d=float(self.ylim_m_box.get())
        except:
            self.error_text["text"] = "Imaginary Min Bound not a number"
            return
        if(a>=b):
            self.error_text["text"] = "Real Bound not valid"
            return
        if(c>=d):
            self.error_text["text"] = "Imaginary Bound not valid"
            return
        area=[a,b,c,d]

        userInput=self.entry_box.get()
        if(ValidateInput(userInput,area)):
            global userFunction
            userFunction=userInput
        else:
            self.error_text["text"] = "Not Valid Function"
            return
            
        Z=RegionLoop(tol,area)
        if(isinstance(Z,str)):
            self.error_text["text"] = Z
            return
        else:
            self.error_text["text"] = ""
        
        Z_=set()
        for n in range(0,len(Z)-1,4):
            x=complex((Z[n]+Z[n+1])/2,(Z[n+2]+Z[n+3])/2);
            x=RoundSigFig(RoundTol(x,tol))
            Z_.add(x)
        
        output=""
        for n in Z_:
            output+=str(n)+"\n"
        
        self.output_area.delete('1.0', tk.END)
        self.output_area.insert(tk.INSERT,output)
        
        Z=list(Z_)
        x=[];
        y=[]
        for z in Z:
            x.append(z.real);
            y.append(z.imag);
        
        fig1 = Figure(figsize=(5,5), dpi=70)
        fig = fig1.add_subplot(111)
        fig.scatter(x,y)

        canvas = FigureCanvasTkAgg(fig1,self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=6, columnspan=4, pady=[0,40])
        
def ValidateInput(userInput,area):
    a=area[0];
    b=area[1];
    c=area[2];
    d=area[3];
    if(userInput=="0"):
        return 0    
    if('_' in userInput):
        return 0
    if('import' in userInput):
        return 0
    x=a+c*1J
    try:
        complex(eval(userInput))
    except:
        return 0
    x=a+d*1J
    try:
        complex(eval(userInput))
    except:
        return 0
    x=b+c*1J
    try:
        complex(eval(userInput))
    except:
        return 0
    x=b+d*1J
    try:
        complex(eval(userInput))
    except:
        return 0
    return 1

def f(x):   #returns "zero" is f(x)=0 else returns f'(x)/f(x)
    try:
        func_x=complex(eval(userFunction))
        if(func_x==0):
            return ("Zero at "+str(x))
        else:
            h=1e-8
            x+=h
            func_xh=complex(eval(userFunction))
            return (func_xh-func_x)/(h*func_x)
    except:
        return ("Error in function at "+str(x))

def RegionLoop(tol,area):
    maxIt=10;
    maxIt_Large=10000;
    Z=[]

    for n in range (0,maxIt_Large):
        region=area[0:4]
        area[0:4]=[]
        for m in range(1,maxIt):
            if(CornorsGood(region)==1):
                num=ComplexRound(ContourIntegral(region,tol)/(pi*1J));
                if(isinstance(num,str)):
                    return num
                if(num.imag==0):
                    break
            region=ScaleRegion(region)
        else:
            return "Scaled region by ~10% and still errored"
        if(num.real>0):
            a=region[0]
            b=region[1]
            c=region[2]
            d=region[3]
            if(abs(a-b)<tol and abs(c-d)<tol):
                Z=Z+region
            else:
                    if(abs(a-b)<abs(c-d)):
                        area=area+[a,b,c,(c+d)/2,a,b,(c+d)/2,d]
                    else:
                        area=area+[a,(a+b)/2,c,d,(a+b)/2,b,c,d]
        if(len(area)==0):
            break
    else:
        return "Looped too many times, tol too small for area"
    return Z

def ComplexRound(x):
    return complex(round(x.real),round(x.imag))

def RoundTol(Z,tol):
    return ComplexRound(Z/tol)*tol

def RoundSigFig(Z):
    s=3
    if(Z.real==0):
        x=0
    else:
        x=round(Z.real,-(math.floor(log10(Z.real).real)-s))
    if(Z.imag==0):
        y=0
    else:
        y=round(Z.imag,-(math.floor(log10(Z.imag).real)-s))
    return complex(x,y) 

def ScaleRegion(area):
    scaleFactor=1.01;

    a=(area[0]+area[1]+(area[0]-area[1])*scaleFactor)/2;
    b=(area[0]+area[1]+(area[1]-area[0])*scaleFactor)/2;
    c=(area[2]+area[3]+(area[2]-area[3])*scaleFactor)/2;
    d=(area[2]+area[3]+(area[3]-area[2])*scaleFactor)/2;

    return [a,b,c,d];

def CornorsGood(area):
    for n in range(0,3):
        if(isinstance(f(area[n]),str)):
            return 0;
    return 1

def ContourIntegral(area,tol):
    a=area[0]
    b=area[1]
    c=area[2]
    d=area[3]

    #new functions of f after parametization
    def f1(t):
        y=f(a+(b-a)*t+c*1J)
        if(isinstance(y,str)):
            return y;
        else:
            return y*(b-a)
    def f2(t):
        y=f(b+1J*((d-c)*t+c));
        if(isinstance(y,str)):
            return y;
        else:
            return y*(d-c)
    def f3(t):
        y=f(b+(a-b)*t+d*1J);
        if(isinstance(y,str)):
            return y;
        else:
            return y*(a-b)
    def f4(t):
        y=f(a+1J*((c-d)*t+d));
        if(isinstance(y,str)):
            return y;
        else:
            return y*(c-d)

    #numerically integrates those functions
    y1=AdaptiveIntegral(f1);
    if(isinstance(y1,str)):
        return y1
    y2=1J*AdaptiveIntegral(f2);
    if(isinstance(y2,str)):
        return y2
    y3=AdaptiveIntegral(f3);
    if(isinstance(y3,str)):
        return y3
    y4=1J*AdaptiveIntegral(f4);
    if(isinstance(y4,str)):
        return y4
    
    return y1+y2+y3+y4;
    

def AdaptiveIntegral(g):
    Nmax=1000
    h_min=1e-10
    h_max=1
    
    tol=0.001;
    h=tol
    t=[0]*(Nmax+1)
    x=[0]*(Nmax+1)

    for n in range(0,Nmax):
        #RK4
        g0=g(t[n])
        if(isinstance(g0,str)):
            return g0;
        g1=g(t[n]+h/4)
        if(isinstance(g1,str)):
            return g1;
        g2=g(t[n]+h/2)
        if(isinstance(g2,str)):
            return g2;
        g3=g(t[n]+3*h/4)
        if(isinstance(g3,str)):
            return g3;
        g4=g(t[n]+h)
        if(isinstance(g4,str)):
            return g4;
        
        y1=x[n]+(g0+4*g2+g4)*h/6;
        y2=x[n]+(g0+4*g1+2*g2+4*g3+g4)*h/12;
        if(abs(y2-y1)>tol): #if RK4 unsuccessful
            h=h*abs(tol/(y2-y1))**0.2
            if(h<h_min):
                h=h_min
            if(h>h_max):
                h=h_max
            g0=g(t[n])
            if(isinstance(g0,str)):
                return g0;
            g2=g(t[n]+h/2)
            if(isinstance(g2,str)):
                return g2;
            g4=g(t[n]+h)
            if(isinstance(g4,str)):
                return g4;
            x[n+1]=x[n]+(g0+4*g2+g4)*h/6
            t[n+1]=t[n]+h
        else:   #if RK4 successful
            x[n+1]=y2
            t[n+1]=t[n]+h
            if(y2-y1==0):
                h=h_max
            else:
                h=h*abs(tol/(y2-y1))**0.2
            if(h<h_min):
                h=h_min
            if(h>h_max):
                h=h_max
        if(t[n+1]>1):
            break
    if(t[n+1]<1):
        print(tol)
        return "Badly behaved function"
    else:
        return x[n]+(1-t[n])/(t[n+1]-t[n])*(x[n+1]-x[n]);

            
root = tk.Tk()
app = Application(master=root)
app.mainloop()
