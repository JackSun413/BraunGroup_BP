import xlwt
import os

SurfaceArea = 0.28274
# lpl_name = str(input("Please enter the folder name:"))
lpl_name = "PLBP_F9ED1_e_1THF_lithiation"

wb = xlwt.Workbook()
ws = wb.add_sheet('the last cycle')


def get_sheet(ssj):
    f12 = f1.readlines()
    lineI_1 = f12[0].split()
    if lineI_1[11] == "Ecell/V" and lineI_1[28] == "Capacity/mA.h":
        i = -1
        x = 0
        y = 0
        lineI = f12[i].split()
        mode = int(lineI[0])
        Ns = int(lineI[6])
        while Ns == 2:
            E_cell = float(lineI[11])
            Specific = float(lineI[28]) / SurfaceArea
            if mode == 3:
                pass
            else:
                ws.write(x, ssj + 3, E_cell)
                ws.write(x, ssj + 2, Specific)
                x += 1

            i += -1
            lineI = f12[i].split()
            Ns = int(lineI[6])
            mode = int(lineI[0])

        while Ns == 1:
            E_cell = float(lineI[11])
            Specific = float(lineI[28]) / SurfaceArea
            if mode == 3:
                pass
            else:
                ws.write(y, ssj + 1, E_cell)
                ws.write(y, ssj, Specific)
                y += 1

            i += -1
            lineI = f12[i].split()
            Ns = int(lineI[6])
            mode = int(lineI[0])

    else:
        # record three important values
        voltage = [o for o, h in enumerate(lineI_1) if h == 'Ecell/V']
        capacity = [o for o, h in enumerate(lineI_1) if h == 'Capacity/mA.h']

        col_heads = ['Ecell/V', 'Capacity/mA.h']
        for cols, ch in zip([voltage, capacity], col_heads):
            assert len(cols) > 0, f'"{ch}" not found in column headers'

        vol = voltage[0]
        capa = capacity[0]
        raw_data = f12[1:]
        capacity, voltage = [], []
        line = f12[number_header_lines].split("\t")
        ecell = float(line[ImZ_col])
        t1me = float(line[vol])
        for line in raw_data:
            each = line.split('\t')
            a = float(each[vol])
            b = float(each[capa]) / SurfaceArea
            if b != 0:
                voltage.append(float(a))
                capacity.append(float(b))


#  D:\科研\data\PLBP_eDEP_6_0_25M_DEE_lithiation\0_25MinDEE_lithiation_13510C_11_GCPL_CA2_IxnQE.txt

wpath = '/Users/sunshijie/Desktop/1/' + lpl_name
files = os.listdir(wpath)

ai = 0
while ai <= 7:
    if not files[ai].find("Cby201053&13510C_03") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(0)
        f1.close()

    if not files[ai].find("Cby201053&13510C_04") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(4)
        f1.close()

    if not files[ai].find("Cby201053&13510C_05") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(8)
        f1.close()

    if not files[ai].find("Cby201053&13510C_06") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(12)
        f1.close()

    if not files[ai].find("Cby201053&13510C_07") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(16)
        f1.close()

    if not files[ai].find("Cby201053&13510C_08") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(20)
        f1.close()

    if not files[ai].find("Cby201053&13510C_09") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(24)
        f1.close()

    if not files[ai].find("Cby201053&13510C_10") == -1:
        f1 = open(wpath + "/" + files[ai], 'r', encoding="latin-1")
        get_sheet(28)
        f1.close()

    ai += 1

wb.save(lpl_name + ".xls")
