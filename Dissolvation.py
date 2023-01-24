import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy.optimize import curve_fit
import cmath
import xlwt

from scipy.optimize import fsolve, curve_fit, least_squares
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d
from scipy.misc import derivative
import scipy.constants as const
import sympy as sp

import pandas as pd
from impedance import preprocessing
from impedance.models.circuits import CustomCircuit
from impedance.visualization import plot_nyquist
import matplotlib.font_manager as font_manager


def read(filename):
    """The purpose of this function is to read the data we need fast.
    But remember, I made two changes to the original data.
    1. The value of current is constant because I read the controlled value of current, not the measured.
    If we want the measured value, you just need to turn h == control/mA to h == I/mA.
    2. The time, even though I record all the time, I didn't adopt it.
    Because if we do the discrete laplace, the default is the time difference is the same all the time.
    So I set the default to be 0.0008.
    If the average time difference is bigger than this, then I will adopt the average number!"""

    # open the file and encoding with latin-1
    with open(filename, 'r', encoding="latin-1") as input_file:
        lines = input_file.readlines()

    header_line = lines[1]

    # MPT data format has variable number of header lines
    number_header_lines = int(header_line.split(":")[1])

    # find the freq and Z columns
    headers = lines[number_header_lines - 1].split('\t')

    # record three important values
    time = [o for o, h in enumerate(headers) if h == 'time/s']
    current = [o for o, h in enumerate(headers) if h == 'control/mA']
    Ecell_2l = [o for o, h in enumerate(headers) if h == '<Ewe>/V']

    col_heads = ['time/s', 'control/mA', '<Ewe>/V']
    for cols, ch in zip([time, current, Ecell_2l], col_heads):
        assert len(cols) > 0, f'"{ch}" not found in column headers'

    freq_col = time[0]
    ReZ_col = current[0]
    ImZ_col = Ecell_2l[0]
    raw_data = lines[number_header_lines + 1:]
    Ecell_2l, current, time = [], [], []
    line = lines[number_header_lines].split("\t")
    ecell = float(line[ImZ_col])
    t1me = float(line[freq_col])
    for line in raw_data:
        each = line.split('\t')
        a = float(each[freq_col]) - t1me
        b = each[ReZ_col]
        c = float(each[ImZ_col]) - ecell
        if b != 0:
            time.append(float(a))
            current.append(float(b))
            Ecell_2l.append(float(c))

    # make some changes to time, make it more concise
    length = len(Ecell_2l)
    tim1 = time[-1]
    delta = tim1 / length
    deviation = (delta - 0.0008) / 0.0008
    if deviation >= 0.1:
        deltaT = delta
    else:
        deltaT = 0.0008
    time1 = length * deltaT

    return Ecell_2l, current, time1, length, deltaT


def function(up, Ecell_2l, current, do=0, gap=1, curve=None, deltat=0.0008):
    """The purpose of this function is to process the read data!
    We normally only have three inputs: the upper limit, original voltage, original current.
    The lower limit(do): If you type nothing, it will start with 0, otherwise just type in whatever number you like.
    Down_sampling(gap): if you want to change the gap, then just type in whatever number you want to change.
    (curve): If You type "curve", then it will curve, otherwise it will output original data.
    Output: New voltage, New current, Frequency"""

    global percentage, Ecell_3l, maximum_deviation
    # we are creating two new list to put the current and voltage into them
    length = (up - do) // gap
    mod = (up - do) % gap
    new_Ecell_2l = []
    new_current = []
    new_time = []
    t1 = (up - do + 1) * deltat
    for i in range(length):
        lpl_E = 0
        lpl_I = 0
        for j in range(gap):
            lpl_E += Ecell_2l[do + i * gap + j]
            lpl_I += current[do + i * gap + j]
        time_0 = i * gap * deltat
        new_time.append(time_0)
        lpl_e = lpl_E / gap
        new_Ecell_2l.append(lpl_e)
        lpl_i = lpl_I / gap
        new_current.append(lpl_i)
    """
    # calculate the last part of the data but I think we can neglect them
    lpl_E = 0
    lpl_I = 0
    for i in range(mod):
        lpl_E += Ecell_2l[up-i]
        lpl_I += current[up-i]
    lpl_e = lpl_E / mod
    new_Ecell_2l.append(lpl_e)
    lpl_i = lpl_I / mod
    new_current.append(lpl_i)
    """
    # curve the function
    if curve == "curve":
        new_ecell = np.array(new_Ecell_2l[10:])
        Ecell_3ll = scipy.signal.savgol_filter(new_ecell, 99, 8)
        Ecell_3l = []
        for i in new_Ecell_2l[:10]:
            Ecell_3l.append(i)
        for i in Ecell_3ll:
            Ecell_3l.append(i)

        # Percentage of deviation
        percentage = []
        for i in range(length):
            a = (Ecell_3l[i - do] - new_Ecell_2l[i - do]) * 100 / Ecell_2l[i - do]
            percentage.append(a)
        percent = get_abs(percentage)
        maximum_deviation = max(percent)

    else:
        # If you do not plan to curve the function using raw data, please use the code below:
        Ecell_3l = new_Ecell_2l
        maximum_deviation = 0

    frequency1 = []
    # frequency calculation
    for i in range(do, up):
        freq = 1 / t1 * (i + 1)
        frequency1.append(freq)
    # return three things, the new_current, the new voltage and the new frequency and the total time
    return Ecell_3l, new_current, frequency1, t1, maximum_deviation


def discrete_laplace_transform(Ecell_3l):
    """The purpose of this function is to do the discrete_laplace_transform"""

    # voltage part, times 2pi/T
    charging = np.array([])
    for i in range(len(Ecell_3l)):
        ssj = Ecell_3l[i] * np.exp(- 2 * np.pi / len(Ecell_3l) * i)
        charging = np.append(charging, [ssj])

    # calculate the Fourier transform
    a = np.fft.fft(charging)
    return a


def expansion(charge, I_current, t1, order):
    """The purpose of this function is to do the Maclaurin expansion"""

    global b

    # get the z(alpha + jw) function
    z_w = []
    for i in range(len(charge)):
        zl = charge[i] / I_current[i]
        z_w.append(zl)

    # get the derivative of the function z_w
    Dydx = []
    Dydx2 = []
    for i in range(len(charge)):
        if i == 0:
            dx = 2 * np.pi / t1 * 1j
            dy = z_w[i] - z_w[i + 1]
            dx2 = dx * dx
            dy2 = z_w[i] - 2 * z_w[i + 1] + z_w[i + 2]
            a = dy / dx
            b = dy2 / dx2
            Dydx.append(a)
            Dydx2.append(b)

        elif i == len(charge) - 1:
            dx = 2 * np.pi / t1 * 1j
            dy = z_w[i - 1] - z_w[i]
            a = dy / dx
            Dydx.append(a)
            Dydx2.append(b)

        else:
            dx = 4 * np.pi / t1 * 1j
            dy = z_w[i - 1] - z_w[i + 1]
            dx2 = dx * dx / 4
            dy2 = z_w[i - 1] - 2 * z_w[i] + z_w[i + 1]
            a = dy / dx
            b = dy2 / dx2
            Dydx.append(a)
            Dydx2.append(b)

    # calculate the z_final(zw) and put them in a list called "z_final"
    z_final = np.array([])

    # different order calculation
    if order == 1:
        for i in range(len(charge)):
            zw = (z_w[i] + Dydx[i] * (-2) * np.pi / t1) * 1000
            z_final = np.append(z_final, [zw])
    elif order == 2:
        for i in range(len(charge)):
            zw = (z_w[i] + Dydx[i] * (-2) * np.pi / t1 + Dydx2 * 2 * np.pi * np.pi / (t1 * t1)) * 1000
            z_final = np.append(z_final, [zw])

    return z_final


def get_abs(x):
    """The purpose of this fuction is to turn all the number in the list/array to positive number"""
    a = abs(np.imag(x))
    b = abs(np.real(x))
    c = []
    for i in range(len(x)):
        d = a[i] * 1j + b[i]
        c.append(d)
    return c


def find_frequency(frequency_list, frequency):
    """The purpose of this function is to find the frequency you want in the original frequency dataset.
    Input: frequency_list means you can put all the number in a list / you can straightly put a number.
    frequency means you need to input the original frequency dataset.
    Output is list anyway: a list with all the numbers in paired frequency_list"""

    order_list = []
    try:
        for i in frequency_list:
            try:
                a = frequency.index(i)
            except:
                D_percent = []
                for j in frequency:
                    lpl = i - j
                    lpl_abs = abs(lpl)
                    D_percent.append(lpl_abs)
                lpl_min = min(D_percent)
                a = D_percent.index(lpl_min)
            order_list.append(a)
    except:
        D_percent = []
        for j in frequency:
            lpl = frequency_list - j
            lpl_abs = abs(lpl)
            D_percent.append(lpl_abs)
        lpl_min = min(D_percent)
        a = D_percent.index(lpl_min)
        order_list.append(a)
    return order_list


def get_sheet_complex(name, num):
    """The purpose of this function is to save the complex number in the list
    the left column is the real number and the right column is the imaginary number"""

    re = np.real(name)
    imag = -np.imag(name)
    for k in range(len(name)):
        ws.write(k, num, re[k])
        ws.write(k, num + 1, imag[k])


def get_sheet(name, num):
    """The purpose of this function is to save the number in the list"""

    for k in range(len(name)):
        ws.write(k, num, name[k])


def show_nyquist(z_final, frequency=None, name=None):
    """The purpose of this function is to get the z_final and turn the imaginary part to positive
    Note: this is different from function get_abs, because some imaginary number might be positive in the first place"""

    real = (np.real(z_final))
    imaginary = -(np.imag(z_final))
    return real, imaginary


def concise(E, I, T, order=1):
    """The purpose of this function is to pack the two functions above and straightly get the z_final from the E and I.
    The other two inputs are the approximate order and the time.
    I strongly recommend you to read all the note for function:
    discrete_laplace_transform and expansion.
    This function will return the z_final directly"""

    # read the text first, to get original data
    E_after = discrete_laplace_transform(E)
    I_after = discrete_laplace_transform(I)

    # Final calculation
    # input: E_after, I_after, Time_total, Order(the order of the simulation)
    z_final = expansion(E_after, I_after, T, order)
    return z_final


def deceitful(number_for_the_first_semicircle, interception):
    """The purpose of this function is to fit better by cutting down the number for the second semicircle.
    We select all the points(number for the first semicircle).
    For the second semicircle, we can select some points. for example,
    if interception = 50, it means for every 50 points, we select one point."""

    Z_new = []
    frequencies_new = []
    testing = 0
    while testing <= number_for_the_first_semicircle:
        Z_new.append(Z[testing])
        frequencies_new.append(frequencies[testing])
        testing += 1
    while testing < len(Z):
        t = testing % interception
        if t == 0:
            Z_new.append(Z[testing])
            frequencies_new.append(frequencies[testing])
            testing += 1
        else:
            testing += 1
    Z_new = np.array(Z_new)
    frequencies_new = np.array(frequencies_new)
    return Z_new, frequencies_new


if __name__ == "__main__":
    # Enter the file name
    name_text = "2DEE_1C_delithiation"

    # Enter the path
    read_text = "/Users/sunshijie/Desktop/Research/Paul Braun(British Petroleum)/DCEIS/dceis data processing/"+ name_text +".mpt"

    # If and only if you type "curve", the system will curve
    choice = 0

    # you can check the notation for function deceitful for detail explanation
    inter = 1

    # guess the number for the circuit, normally the relative value matter, not the absolute value
    R0 = 1.70e+02
    R1 = 6.27e+01
    CPE_1_0 = 4.18e-03
    CPE_1_1 = 1
    R2 = 5.79e+02
    CPE_2_0 = 1.38e-05
    CPE_2_1 = 7.89e-01

    wb = xlwt.Workbook()
    ws = wb.add_sheet('Sheet 1')

    Ei, Ii, T, N, deltaT = read(read_text)
    E, I, F, t1, maximum_deviation = function(N, Ei, Ii, deltat=deltaT, curve=choice)

    # Do the Laplace transform and get the original data
    z_final = concise(E, I, t1, order=1)

    # save all the origin data for z_final
    Num = int(N / 2)
    get_sheet(F[1:Num], 6)
    get_sheet_complex(z_final[1: Num], 7)

    # find the maximal and minimal for all the data
    Mini = min(-np.imag(z_final[10:400]))
    z_final1 = list(-np.imag(z_final[0:400]))
    mini = z_final1.index(Mini)
    mi = F[mini]
    Maxi = max(-np.imag(z_final[400:2500]))
    z_final2 = list(-np.imag(z_final[0:2500]))
    maxi = z_final2.index(Maxi)
    ma = F[maxi]

    # you can put all the data you want to specially annotate in the list
    typical = [0.1, 1, 10, 100, ma]

    # you can change the frequency value, this will determine how many data you want to draw and fit
    frequency_you_want_to_fit = 100

    Z = []
    frequencies = []

    # this is to prevent you enter something other than the typical number in the list typical
    try:
        find_100 = typical.index(frequency_you_want_to_fit)
    except:
        typical_transition = []
        for i in typical:
            typical_transition.append(i)
        typical_transition.append(frequency_you_want_to_fit)
        typical = typical_transition
        find_100 = typical.index(frequency_you_want_to_fit)

    order_list = find_frequency(typical, F)
    order_100 = order_list[find_100]
    for i in range(0, order_100):
        Z.append(z_final[i])
        frequencies.append(F[i])
    Z = np.array(Z)  # Turn the impedance list into an array
    frequencies = np.array(frequencies)

    get_sheet_complex(Z, 1)  # zw
    get_sheet(frequencies, 0)  # frequency

    # C# is a capacitance, Wo# is a warburg element, etc. There are more in the doc
    # https://impedancepy.readthedocs.io/en/latest/circuit-elements.html
    circuit = 'R0-p(R1,CPE_1)-p(R2,CPE_2)'  # Choose your model R# is a resistance, - means series, p(,) means parallel

    # Set your initial guess for each element
    initial_guess = [R0, R1, CPE_1_0, CPE_1_1, R2, CPE_2_0, CPE_2_1]

    # Create the circuit with those guess parameters
    circuit = CustomCircuit(circuit, initial_guess=initial_guess)

    # Note!!! Here we change the original data for fitting
    # If you want to use the original data, please delete next line
    Z, frequencies = deceitful(100, inter)

    # Fit the circuit to your data
    circuit.fit(frequencies, Z)

    # This outputs your guess values and the fitting results with their units
    print(circuit)

    # save the fitting data: zw
    Z_fit = circuit.predict(frequencies)
    get_sheet_complex(Z_fit, 3)
    wb.save(name_text + ".xls")

    # Below are the drawing part!!!
    # Below are the drawing part!!!
    # Below are the drawing part!!! Please don't change anything!!!
    fig, ax = plt.subplots()

    # annotate all the number you want to specially note
    note = 0
    note_max = len(typical)
    while note < note_max:
        i = order_list[note]
        j = typical[note]
        special_y = -np.imag(z_final[i])
        special_x = np.real(z_final[i])
        plt.annotate("%d Hz" % j, (special_x, special_y), xytext=(special_x, Maxi * 1.1), arrowprops=dict(arrowstyle="->"), fontname="Arial",
                     fontsize=14)
        note += 1

    # this set the basic property for the plot
    font = font_manager.FontProperties(family='Arial', style='normal', size=18)
    fig.set_figheight(8)
    fig.set_figwidth(10)

    # this plots your actual data. lw controls the thickness of lines
    plot_nyquist(ax, Z, color='blue', fmt='o', label='LT Impedance Data', lw=2.5)

    # this plots your fit data
    plot_nyquist(ax, Z_fit, color='orange', fmt='-', label='Fitting', lw=2.5)

    # these three lines make the figure looks square
    length_limit = np.real(Z[0]) * 1.1
    ax.set_ylim(0, length_limit)  # Set the y-axis limits
    ax.set_xlim(0, length_limit)  # Set the x-axis limits

    # these two lines label two axis
    ax.set_xlabel("Z'(ω) [Ω]", fontname="Arial", fontsize=19)
    ax.set_ylabel('-Z"(ω) [Ω]', fontname="Arial", fontsize=19)

    # plt.setp(ax.get_xticklabels(), fontsize=18, fontname="Arial") #Set the font of your ticks if they are numbers
    # plt.setp(ax.get_yticklabels(), fontsize=18, fontname="Arial") #Set the font of your ticks if they are numbers

    ax.tick_params(axis="y", direction="in")  # Direction of ticks
    ax.tick_params(axis="x", direction="in")  # Direction of ticks
    # ax.set_yticklabels([],fontname='Arial',fontsize=25) #You can choose which values are shown

    plt.legend(prop=font, loc='upper right')  # Control position of legend and the font which is set earlier
    ax.set_axisbelow(True)  # Forces data to above any figure items
    ax.grid(linestyle='dotted')  # Sets a grid and choose the linestyle
    plt.show()

