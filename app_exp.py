import tkinter as tk
from tkinter import ttk
import tkinter.font as tkFont
from tkinter import filedialog as fd
import numpy as np
import matplotlib
import scipy.stats
from scipy.signal import find_peaks
from tkinter.messagebox import showerror, showwarning, showinfo
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
from matplotlib.collections import LineCollection



def find_zero(array, array1, direction, nmax, tol1, tol2):
    for ind in range(1,nmax):
        if direction=='left':
            ind *= -1
        v = array1[ind]
        if abs(v) < tol1 and abs(array[ind])<abs(tol2):
            return abs(ind)
    return 0
m = 1e-6
class Signal:
    def __init__(self, fname):
        self.fname = fname
        self.impath = None
        self.vrods = None
        self.vmesh = None
        self.x = None
        self.y = None
        self.y1 = None
        self.ca_int = []
        self.bl_int = []
        self.bl_dx = 80000
        self.bl_center = -25000
        self.y_smooth = None
        self.peak_ind = []
        self.peak_time = []
        self.peak_amp = []
        self.peak_start = []
        self.peak_end = []
        self.bl_check = []
        
    
    def get_voltages(self): #if needed
        self.vrods = self.fname[29:32]
        self.vmesh = self.fname[38:41]

    def get_data(self):
        data = np.loadtxt(self.fname)
        self.x, self.y = data[:,0], data[:,2]
        
    def smooth_y(self, w_length=51, polyorder=3):
        self.y_smooth = savgol_filter(self.y, w_length, polyorder) #x, window_length (must be odd), polyorder

    def get_peak_data(self, npeak=1, height=0.8, width=20, rel_height=0.3, threshold=3*10e-5):
        #If the peaks are not well defined, probably change the threshold
        self.peak_ind = find_peaks(-self.y_smooth, height=height)[0][0:npeak]
        self.peak_time = self.x[self.peak_ind]
        self.peak_amp = self.y_smooth[self.peak_ind]

        self.y1 = (np.diff(self.y_smooth))
        self.peak_start = np.zeros(len(self.peak_ind), dtype=int)
        self.peak_end = np.zeros(len(self.peak_ind), dtype=int)
        for i, p in enumerate(self.peak_ind):
            self.peak_start[i] = p - find_zero(self.y1[p-500:p-10], 'left', 300, threshold)
            self.peak_end[i] = p + find_zero(abs(self.y1[p+10:]), 'right', 300, threshold)

    def get_integrals(self, shift=1, trapz=True):
        dx = self.x[1] - self.x[0] #spacing between 2 points in time
        self.ca_int = np.zeros(len(self.peak_ind))
        self.bl_int = np.zeros(len(self.peak_ind))
        self.bl_check = np.zeros(len(self.peak_ind))
        for i, p in enumerate(self.peak_ind):
            ps, pe = self.peak_start[i], self.peak_end[i]
            bl_dx = pe-ps #Window of integration
            if trapz:
                
                self.ca_int[i] = np.trapz(-self.y[ps:pe], self.x[ps:pe])
                self.bl_int[i] = np.trapz(-self.y[-bl_dx-shift:-shift], self.x[-bl_dx-shift:-shift])
            else:
                #Rectangle method
                self.ca_int[i] = np.sum(-self.y[ps:pe])*dx
                self.bl_int[i] = np.sum(-self.y[-bl_dx-shift:-shift])*dx
                self.bl_check[i] =(np.sum(-self.y[ps:pe]) - np.sum(-self.y[-bl_dx-shift:-shift]))*dx
                # print(i, ps, pe, self.ca_int[i])

            
    def norm_x(self): #in ms instead of s
        self.x /= m
        self.y = self.y[self.x>0]
        self.x = self.x[self.x>0]
        

class App(tk.Tk):
    def exit(self):
        self.quit()
        self.destroy()


    def __init__(self):
        self.y_min = 1
        self.is_normed = False
        super().__init__()
        self.protocol('WM_DELETE_WINDOW', self.exit)
        default_font = tkFont.nametofont('TkDefaultFont')
        default_font.configure(size=13, family='Comic Sans MS')

        # self.geometry("240x100")
        self.title('Astrochemistry data analyzer')
        # self.resizable(0, 0)

        self.left_frame = ttk.Frame(self, padding=1, borderwidth = 5)
        self.plot_frame = ttk.Frame(self, padding=1)
        self.right_frame = ttk.Frame(self, padding=1)
        self.int_frame = ttk.Frame(self, padding=1)

        self.create_left_frame(self.left_frame)
        self.create_plot_frame(self.plot_frame)
        self.create_right_frame(self.right_frame)
        
        self.left_frame.pack(side=tk.LEFT, anchor=tk.N)
        self.plot_frame.pack(side=tk.LEFT, anchor=tk.N)


    def create_left_frame(self, left_frame):
        #Choose exp or sim
        type_frame = ttk.Frame(left_frame)
        
        file_type_label = ttk.Label(type_frame, text="Type of data :")
        self.file_cbox= ttk.Combobox(type_frame)
        self.file_cbox['values'] = ['Experiment']#['Experiment', 'Simulation']
        self.file_cbox['state'] = 'readonly'
        self.file_cbox.current(0)

        file_button = ttk.Button(left_frame, text='Open a File', command=self.select_file)
        self.norm_button = ttk.Button(left_frame, text='Normalize')
        self.unnorm_button = ttk.Button(left_frame, text='Undo normalization')
        self.bl_button = ttk.Button(left_frame, text='Adjust baseline')
        # self.clean_btn = ttk.Button(left_frame, text='Clean window', command=self.clean)

        
        self.cal1_frame = ttk.Frame(left_frame)
        self.calib_label = ttk.Label(self.cal1_frame, text='\n\nCalibration for m/z axis')
        self.masslabel = ttk.Label(self.cal1_frame, text='known m/z')
        self.toflabel = ttk.Label(self.cal1_frame, text='known TOF')
        self.mass1_entry = ttk.Entry(self.cal1_frame)
        self.mass2_entry = ttk.Entry(self.cal1_frame)
        self.tof1_entry = ttk.Entry(self.cal1_frame)
        self.tof2_entry = ttk.Entry(self.cal1_frame)
        self.calib_btn = ttk.Button(self.cal1_frame, text='Calibration', command=self.calib)

        self.calib_label.grid(row=0, columnspan=2)
        self.masslabel.grid(row=1, column=0)
        self.toflabel.grid(row=1, column=1)
        self.mass1_entry.grid(row=2, column=0)
        self.mass2_entry.grid(row=3, column=0)
        self.tof1_entry.grid(row=2, column=1)
        self.tof2_entry.grid(row=3, column=1)
        self.calib_btn.grid(row=4, columnspan=2)

        type_frame.pack()
        file_type_label.pack(side=tk.LEFT)
        self.file_cbox.pack(side=tk.LEFT)
        file_button.pack()
        # self.clean_btn.pack()
        # self.cal1_frame.pack()

    def calib(self):
        def mass_to_tof1(x):
            return known_tof1*np.sqrt(x/known_mass1)
        def tof_to_mz1(x):
            k = known_mass1 / (known_tof1**2)
            return k * x**2
        try:
            self.secax.remove()
        except Exception as e:
            pass
        if self.mass1_entry.get() != '' and self.tof1_entry.get() != '':
            known_mass1 = float(self.mass1_entry.get())
            known_tof1 = float(self.tof1_entry.get())
            
            if self.mass2_entry.get() != '' and self.tof2_entry.get() != '':
                known_mass2 = float(self.mass2_entry.get())
                known_tof2 = float(self.tof2_entry.get())
                known_tof = np.array([known_tof1, known_tof2])
                known_mass = np.array([known_mass1, known_mass2])
                coefficients = np.polyfit(known_tof**2, known_mass, 1)
                def tof_to_mz2(x):
                    return coefficients[0] * x**2 + coefficients[1]
                coefficients_m = np.polyfit(np.sqrt(known_mass), known_tof, 1)
                def mass_to_tof2(x):
                    return coefficients_m[0] * np.sqrt(x) + coefficients_m[1]
                self.secax = self.main_ax.secondary_xaxis('top', functions=(tof_to_mz2,mass_to_tof2))
                self.secax.set_xlabel('m/z')
                self.main_figure_canvas.draw()
            else:
                self.secax = self.main_ax.secondary_xaxis('top', functions=(tof_to_mz1,mass_to_tof1))
                self.secax.set_xlabel('m/z')
                self.main_figure_canvas.draw()
        else:
            print('no cond')


    def clean(self):
        Frames = [self.plot_frame, self.left_frame, self.right_frame, self.int_frame]
        for frame in Frames:
            for widgets in frame.winfo_children():
                widgets.destroy()
        self.create_left_frame(self.left_frame)
        self.create_plot_frame(self.plot_frame)
        self.create_right_frame(self.right_frame)
        self.right_frame.pack_forget()
        self.int_frame.pack_forget()

    def create_int_frame(self, int_frame):
        pass

    def create_plot_frame(self, plot_frame):
        pass


    def create_right_frame(self, right_frame):

        h_label = ttk.Label(right_frame, text="Min. height (abs)")
        self.h = tk.DoubleVar(right_frame, value=0.1)
        h_entry = ttk.Entry(right_frame, textvariable=self.h)
        
        dist_label = ttk.Label(right_frame, text="Distance between peaks (micros)")
        self.dist = tk.DoubleVar(right_frame, value=0.12)
        dist_entry = ttk.Entry(right_frame, textvariable=self.dist)
        

        self.autopeak_button = ttk.Button(right_frame, text='Autofind peaks')
        self.autoint_button = ttk.Button(right_frame, text='Auto integration')

        int_label = ttk.Label(right_frame, text="Integrations limit:")
        x1_int_label = ttk.Label(right_frame, text="t1 (micros)")
        self.x1_int = tk.DoubleVar(right_frame, value=0)
        x1_int_entry = ttk.Entry(right_frame, textvariable=self.x1_int)
        
        x2_int_label = ttk.Label(right_frame, text="t2 (micros)")
        self.x2_int = tk.DoubleVar(right_frame, value=1)
        x2_int_entry = ttk.Entry(right_frame, textvariable=self.x2_int)

        self.man_int_button = ttk.Button(right_frame, text='Manual integration')
        self.man_int_value = tk.DoubleVar(right_frame)
        self.man_int_label = ttk.Label(right_frame)

        self.thr_der_label = ttk.Label(right_frame, text="Threshold on the derivative")
        self.thr_der = tk.DoubleVar(right_frame, value=3e-4)
        self.thr_der_entry = ttk.Entry(right_frame, textvariable=self.thr_der)
        
        self.thr_peak_label = ttk.Label(right_frame, text="1/Threshold on the signal (fraction of the peak)")
        self.thr_peak = tk.DoubleVar(right_frame, value=10)
        self.thr_peak_entry = ttk.Entry(right_frame, textvariable=self.thr_peak)
        
        h_label.pack()
        h_entry.pack()
        dist_label.pack()
        dist_entry.pack()
        x1_int_label.pack()
        x1_int_entry.pack()
        x2_int_label.pack()
        x2_int_entry.pack()
        self.man_int_button.pack()


    def select_file(self):
        if(self.file_cbox.get() == ''):
            showerror(title='Error', message='Please choose between experiment and simulation')
        else:
            filetypes = (
                ('text files', '*.txt'),
                ('All files', '*.*')
            )

            filename = fd.askopenfilename(
                title='Open a file',
                initialdir='./',
                filetypes=filetypes)
            try:
                self.clean()
                self.s = self.load_exp_data(filename)
                numberOfScreenUnits = 200
                self.success_label = ttk.Label(self.left_frame, text=(filename, 'loaded'), font=['Comic Sans MS', 8], wraplength=numberOfScreenUnits)
                self.success_label.pack(side=tk.TOP, anchor=tk.W)

                self.main_fig = plt.figure()
                self.main_figure_canvas = FigureCanvasTkAgg(self.main_fig, self.plot_frame)
                self.main_ax = self.main_fig.add_subplot()
                self.plot_exp_data(self.s, self.main_ax, self.main_figure_canvas)

            except Exception as e:
                showerror(title='Error', message=e)


    def load_exp_data(self, filename):
        s = Signal(filename)
        s.get_data()
        s.norm_x()
        s.smooth_y()
        return s


    def find_peaks(self, s, h, npeak, width=30):
        s.get_peak_data(npeak=npeak.get(), height=h.get(), width=width)
        for i, p in enumerate(s.peak_ind):
            self.main_ax.vlines(s.x[p], 0, min(s.y), color='red', ls='dashed', alpha=0.5)
            self.main_ax.vlines(s.x[s.peak_start[i]], 0, min(s.y), color='pink', ls='dashed')
            self.main_ax.vlines(s.x[s.peak_end[i]], 0, min(s.y), color='pink', ls='dashed')
            self.main_figure_canvas.draw()


    def plot_exp_data(self, s, ax, fig_canvas):
        self.main_line = ax.plot(s.x, s.y)
        ax.set_xlabel('TOF (micros)')
        ax.set_ylabel('Signal (V)')
        toolbar_frame = tk.Frame(self.plot_frame)
        toolbar = NavigationToolbar2Tk(fig_canvas, toolbar_frame) 

        self.norm_button.configure(command=lambda:self.norm(s, ax, fig_canvas))
        self.unnorm_button.configure(command=lambda:self.unnorm(s, ax, fig_canvas))
        self.bl_button.configure(command=lambda:self.adjust_baseline(s, ax, fig_canvas))
        self.autopeak_button.configure( command=lambda: self.autofind_peaks(s, ax, fig_canvas))
        self.man_int_button.configure(command=lambda: self.manual_int(s))
        self.autoint_button.configure(command=lambda: self.auto_int(s))
        self.norm_button.pack()
        self.unnorm_button.pack()
        self.bl_button.pack()

        self.thr_der_label.pack()
        self.thr_der_entry.pack()
        self.thr_peak_label.pack()
        self.thr_peak_entry.pack()
        self.autopeak_button.pack()
        self.autoint_button.pack()
        fig_canvas.get_tk_widget().pack()
        toolbar_frame.pack() 
        self.right_frame.pack(side=tk.LEFT, anchor=tk.N)
        self.int_frame.pack(side=tk.LEFT, anchor=tk.N)
        self.cal1_frame.pack()


    def manual_int(self, s):
        dx = s.x[1]-s.x[0]
        x1, x2 = int(self.x1_int.get()//dx), int(self.x2_int.get()//dx)
        integral = np.trapz(s.y[x1:x2], s.x[x1:x2])
        self.man_int_label.configure(text= '{:.5g}'.format(integral))
        self.man_int_label.pack()

    def auto_int(self, s):
        for widget in self.int_frame.winfo_children():
            widget.destroy()
        dx = s.x[1]-s.x[0]
        subframe = tk.Frame(self.int_frame)
        l1 = tk.Label(subframe, text='Nb of the peak', bd=2)
        l2 = tk.Label(subframe, text='TOF (micros)', bd=2)
        l3 = tk.Label(subframe, text='Start (micros)', bd=2)
        l4 = tk.Label(subframe, text='End (micros)', bd=2)
        l5 = tk.Label(subframe, text='Integral', bd=2)
        l1.grid(row=0, column=0)
        l2.grid(row=0, column=1)
        l3.grid(row=0, column=2)
        l4.grid(row=0, column=3)
        l5.grid(row=0, column=4)
        for i in range(len(s.peak_time)):
            
            nb = tk.Label(subframe)
            nb['text'] = i
            nb.grid(row=i+1, column=0)

            tof = tk.Label(subframe)
            tof['text'] = '{:.5g}'.format(s.peak_time[i])
            tof.grid(row=i+1, column=1)

            ps = tk.Label(subframe)
            ps['text'] = '{:.5g}'.format(s.peak_start[i]*dx)
            ps.grid(row=i+1, column=2)

            pe = tk.Label(subframe)
            pe['text'] = '{:.5g}'.format(s.peak_end[i]*dx)
            pe.grid(row=i+1, column=3)


            x1, x2 = int(s.peak_start[i]), int(s.peak_end[i])
            integral = np.trapz(s.y[x1:x2], s.x[x1:x2])
            integ = tk.Label(subframe)
            integ['text'] = '{:.6g}'.format(integral)
            integ.grid(row=i+1, column=4)
        subframe.grid()


    def norm(self, s, ax, fig_canvas):
        if not self.is_normed:
            self.is_normed = True
            self.y_min = np.min(s.y)
            s.y /= -self.y_min
            # ax.clear()
            line = self.main_line.pop(0)
            line.remove()
            self.main_line =ax.plot(s.x, s.y, color='blue')
            fig_canvas.draw()


    def unnorm(self, s, ax, fig_canvas):
        if self.is_normed:
            self.is_normed = False
            s.y *= -self.y_min
            # ax.clear()
            # del self.main_line
            line = self.main_line.pop(0)
            line.remove()
            self.main_line = ax.plot(s.x, s.y, color='blue')
            fig_canvas.draw()

    def adjust_baseline(self, s, ax, fig_canvas):
        dx = s.x[1]-s.x[0]
        s.y -= scipy.stats.mode(s.y[int(32//dx):int(37//dx)])[0]
        # ax.clear()
        # ax.lines.pop(0)
        line = self.main_line.pop(0)
        line.remove()
        ax.vlines([32,37], 0,-1, ls='dotted', color='black', alpha=0.3)
        self.main_line = ax.plot(s.x, s.y, color='blue')
        fig_canvas.draw()


    def autofind_peaks(self, s, ax, fig_canvas):
        dx = s.x[1]-s.x[0]
        y_copy = np.copy(s.y)
        y_copy /= -np.min(y_copy)
        if int(self.dist.get()//dx) != 0:
            P =  find_peaks(-y_copy, height=self.h.get(), distance=int(self.dist.get()//dx))[0]
        else:
            P =  find_peaks(-y_copy, height=self.h.get(), distance=300)[0]
        P = P[(s.x[P]>1) & (s.x[P]<30)]
        threshold = float(self.thr_der.get())
        y1 = np.diff(y_copy)
        #def find_zero(array, array1, direction, nmax, tol1, peak_h):
        PS = [p -10 - find_zero(y_copy[p-500:p-10], y1[p-500:p-10], 'left', 400, threshold, y_copy[p]/float(self.thr_peak.get())) for p in P]
        PE = [p + 10 + find_zero(y_copy[p+10:], y1[p+10:], 'right', 400, threshold, y_copy[p]/float(self.thr_peak.get())) for p in P]
        s.peak_start = PS
        s.peak_end = PE
        s.peak_ind = P
        s.peak_time = s.x[P]
        try:
            for child in ax.get_children():
                if isinstance(child, LineCollection):
                    child.remove()
        except:
            pass
        ax.vlines(s.x[P], 0,-1,color='red', ls='dashed', alpha=0.5)
        ax.vlines(s.x[PE], 0,-1,color='orange', ls='dashed', alpha=0.5)
        ax.vlines(s.x[PS], 0,-1,color='orange', ls='dashed', alpha=0.5)
        fig_canvas.draw()


if __name__ == "__main__":
    app = App()
    app.mainloop()