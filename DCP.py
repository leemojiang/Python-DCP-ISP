import numpy as np
import xml.etree.ElementTree as ET
import math
# from lxml import etree as ET
from xml.dom import minidom


def indent(elem, level=0):
    # print(elem)
    # https://blog.csdn.net/hiccupfrost/article/details/107594835
    i = "\n" + level*"\t"
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "\t"
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def XYZ2xyY(XYZ):
    '''
        XYZ [X,Y,Z] Shape in [3,]
        xyY [x,y,Y]

        Ref: http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_xyY.html

    '''
    X,Y,Z = XYZ

    x = X / (X+Y+Z)
    y = Y / (X+Y+Z)
    
    return [x,y,Y]

def xy2Temp(x,y):
    '''
        Ref: https://www.waveformlighting.com/tech/calculate-color-temperature-cct-from-cie-1931-xy-coordinates
    '''

    n = (x-0.3320)/(0.1858-y);
    CCT = 437*n**3 + 3601*n**2 + 6861*n + 5517
    return CCT

def xyCoordToTemperature(white_xy):
    '''
        Ref: https://github.com/Beep6581/RawTherapee/blob/cc0e941652b7d0057411ad04299a21203fc3e2bb/rtengine/dcp.cc#L1857

    '''


    class Ruvt:
        def __init__(self, r, u, v, t):
            self.r = r
            self.u = u
            self.v = v
            self.t = t

    temp_table = [
        Ruvt(0, 0.18006, 0.26352, -0.24341),
        Ruvt(10, 0.18066, 0.26589, -0.25479),
        Ruvt(20, 0.18133, 0.26846, -0.26876),
        Ruvt(30, 0.18208, 0.27119, -0.28539),
        Ruvt(40, 0.18293, 0.27407, -0.30470),
        Ruvt(50, 0.18388, 0.27709, -0.32675),
        Ruvt(60, 0.18494, 0.28021, -0.35156),
        Ruvt(70, 0.18611, 0.28342, -0.37915),
        Ruvt(80, 0.18740, 0.28668, -0.40955),
        Ruvt(90, 0.18880, 0.28997, -0.44278),
        Ruvt(100, 0.19032, 0.29326, -0.47888),
        Ruvt(125, 0.19462, 0.30141, -0.58204),
        Ruvt(150, 0.19962, 0.30921, -0.70471),
        Ruvt(175, 0.20525, 0.31647, -0.84901),
        Ruvt(200, 0.21142, 0.32312, -1.0182),
        Ruvt(225, 0.21807, 0.32909, -1.2168),
        Ruvt(250, 0.22511, 0.33439, -1.4512),
        Ruvt(275, 0.23247, 0.33904, -1.7298),
        Ruvt(300, 0.24010, 0.34308, -2.0637),
        Ruvt(325, 0.24702, 0.34655, -2.4681),
        Ruvt(350, 0.25591, 0.34951, -2.9641),
        Ruvt(375, 0.26400, 0.35200, -3.5814),
        Ruvt(400, 0.27218, 0.35407, -4.3633),
        Ruvt(425, 0.28039, 0.35577, -5.3762),
        Ruvt(450, 0.28863, 0.35714, -6.7262),
        Ruvt(475, 0.29685, 0.35823, -8.5955),
        Ruvt(500, 0.30505, 0.35907, -11.324),
        Ruvt(525, 0.31320, 0.35968, -15.628),
        Ruvt(550, 0.32129, 0.36011, -23.325),
        Ruvt(575, 0.32931, 0.36038, -40.770),
        Ruvt(600, 0.33724, 0.36051, -116.45)
    ]

    res = 0

    u = 2.0 * white_xy[0] / (1.5 - white_xy[0] + 6.0 * white_xy[1])
    v = 3.0 * white_xy[1] / (1.5 - white_xy[0] + 6.0 * white_xy[1])

    last_dt = 0.0

    for index in range(1, 31):
        du = 1.0
        dv = temp_table[index].t
        length = math.sqrt(1.0 + dv * dv)
        du /= length
        dv /= length

        uu = u - temp_table[index].u
        vv = v - temp_table[index].v

        dt = -uu * dv + vv * du

        if dt <= 0.0 or index == 30:
            if dt > 0.0:
                dt = 0.0

            dt = -dt

            f = 0.0 if index == 1 else dt / (last_dt + dt)

            res = 1.0e6 / (temp_table[index - 1].r * f + temp_table[index].r * (1.0 - f))
            break

        last_dt = dt

    return res

def interpolateMatrix(temp1,temp2,wbtemp,MA1,MA2):
    '''
        From DNG specification

        Ref: https://github.com/Beep6581/RawTherapee/blob/cc0e941652b7d0057411ad04299a21203fc3e2bb/rtengine/dcp.cc#L1857

    '''

    if wbtemp <= temp1:
        mix = 1.0
    elif wbtemp >= temp2:
        mix = 0.0
    else:
        invT = 1.0 / wbtemp
        mix = (invT - (1.0 / temp2)) / ((1.0 / temp1) - (1.0 / temp2))
    
    
    mixM = mix * MA1 + (1-mix)*MA2
    return mixM

class DCPParser:
    # # Usage
    # parser = DCPMatrixParser('.\Sony ILCE-7RM2 Adobe Standard.xml')
    # parser.parse_matrices()
    # parser.matrices['ColorMatrix1'] = np.zeros((3,3))
    # # Modify matrices if needed
    # parser.save_matrices('./output.xml',tags_to_exclude=['HueSatDeltas1', 'HueSatDeltas2','LookTable'])

    def __init__(self, file_path):
        self.file_path = file_path
        self.matrices = {}
        self.tag_pairs = {}
        self.luts = {}

        self.tree = ET.parse(self.file_path)

    def modify_tag_pairs(self, pairs = {"Copyright": "Not Adobe"}):
        for k,v in pairs.items():
            self.tag_pairs[k] = v
    
    def parse_matrices(self):
        # tree = ET.parse(self.file_path)
        tree =self.tree
        root = tree.getroot()
        for matrix_elem in root.findall('*'):
            if 'Matrix' in matrix_elem.tag:
                matrix_name = matrix_elem.tag
                rows = int(matrix_elem.attrib['Rows'])
                cols = int(matrix_elem.attrib['Cols'])
                matrix = np.zeros((rows, cols))
                for element in matrix_elem.findall('Element'):
                    row = int(element.attrib['Row'])
                    col = int(element.attrib['Col'])
                    value = float(element.text)
                    matrix[row][col] = value
                self.matrices[matrix_name] = matrix

    def parse_LUTs(self):
        # LookTable HueSatDeltas2 HueSatDeltas1
        tree = self.tree
        root = tree.getroot()

        for lut_elem in root.findall('*'):
            if 'LookTable' == lut_elem.tag or 'HueSatDeltas' in lut_elem.tag :
                lut_name = lut_elem.tag
                
                hue_div = int(lut_elem.attrib['hueDivisions']) # ex hueDivisions=90 0,89
                sat_div = int(lut_elem.attrib['satDivisions'])
                val_div = int(lut_elem.attrib['valDivisions'])
                tensor = np.zeros((hue_div, sat_div, val_div,3)) #[h index, s index, v index, hsv] hsv hue shift sv scale

                print(f"{lut_name} dim: {tensor.shape} ")
                
                for element in lut_elem.findall('Element'):
                    h = int(element.attrib['HueDiv'])
                    s = int(element.attrib['SatDiv'])
                    v = int(element.attrib['ValDiv'])

                    h_ = float(element.attrib['HueShift'])
                    s_ = float(element.attrib['SatScale'])
                    v_ = float(element.attrib['ValScale'])   

                    tensor[h,s,v,:] = [h_,s_,v_]
                
                self.luts[lut_name] = tensor


    def save_modified(self, output_file_path,rename_profileName=None ,tags_to_exclude=[]):
        tree = ET.parse(self.file_path)
        root = tree.getroot()

        for k,v in self.tag_pairs.items():
            for element in root.findall(k): # root.find
                previous = element.text
                element.text = v
                print(f"Replace \"{k}\" to {v}")

        if rename_profileName is not None:
            for element in root.findall("ProfileName"): # root.find
                element.text = rename_profileName
                print(f"Rename \"ProfileName\" to {rename_profileName}")

        # Rewrite matrix 
        for matrix_name, matrix in self.matrices.items():
            if matrix_name in tags_to_exclude:
                continue

            new_matrix_elem = ET.Element(matrix_name)
            new_matrix_elem.attrib['Rows'] = str(matrix.shape[0])
            new_matrix_elem.attrib['Cols'] = str(matrix.shape[1])

            #Reverse order to result in xml
            for row in reversed(range(matrix.shape[0])):
                for col in reversed(range(matrix.shape[1])): 
                    element = ET.Element('Element')
                    element.attrib['Row'] = str(row)
                    element.attrib['Col'] = str(col)
                    element.text = str(matrix[row][col])
                    new_matrix_elem.append(element)

             # Remove existing matrix elements
            for elem in root.findall(matrix_name):
                root.remove(elem)
            
            root.append(new_matrix_elem)

        #  # Write Modified Matrix
        # for matrix_elem in root.findall('*'):
        #     if matrix_elem.tag in self.matrices: #keys
        #         matrix = self.matrices[matrix_elem.tag]
        #         for element in matrix_elem.findall('Element'):
        #             row = int(element.attrib['Row'])
        #             col = int(element.attrib['Col'])
        #             element.text = str(matrix[row][col])

        # Write new LUTs and remove pervious one
        for lut_name, tensor in self.luts.items():
            if lut_name in tags_to_exclude:
                continue

            new_lut_elem = ET.Element(lut_name)
            new_lut_elem.attrib['hueDivisions'] = str(tensor.shape[0]) # hue_div
            new_lut_elem.attrib['satDivisions'] = str(tensor.shape[1]) # sat_div
            new_lut_elem.attrib['valDivisions'] = str(tensor.shape[2]) # val_div

        
            for h in range(tensor.shape[0]):
                for s in range(tensor.shape[1]):
                    for v in range(tensor.shape[2]):
                        element = ET.Element('Element')
                        element.attrib['HueDiv'] = str(h)
                        element.attrib['SatDiv'] = str(s)
                        element.attrib['ValDiv'] = str(v)
                        element.attrib['HueShift'] = str(tensor[h, s, v, 0])
                        element.attrib['SatScale'] = str(tensor[h, s, v, 1])
                        element.attrib['ValScale'] = str(tensor[h, s, v, 2])
                        new_lut_elem.append(element)
                 
            # Remove existing LUT elements
            for lut_elem in root.findall('*'):
                if lut_elem.tag == lut_name:
                    root.remove(lut_elem)
            root.append(new_lut_elem)   

        for tag in tags_to_exclude:
            elements = root.findall(tag)
            for element in elements:
                root.remove(element)
                print(f"Remove {element}")

        indent(root)
        tree.write(output_file_path, encoding='utf-8', xml_declaration=True)
        # xml_str = minidom.parseString(ET.tostring(root)).toprettyxml()
        # xml_write = minidom.parseString(ET.tostring(root))
        # with open(output_file_path, 'w') as fh:
        #     xml_write.writexml(fh, indent='\t', newl='\n', encoding='utf-8')

import subprocess
import os

class DCPHandler:
    # Usage 
    # file_path = '.\Sony ILCE-7RM2 Adobe Standard.dcp'

    # dcp = DCPHandler(file_path)
    # dcp.convert2xml()
    # dcp.parse()
    # dcp.parser.parse_matrices()

    # dcp.save_modified("./data/ABC.xml",dcp_outpath='./data',rename_profileName="Test Matrix Only",tags_to_exclude=['HueSatDeltas1','HueSatDeltas2','LookTable'])


    def __init__(self, file_path , dcpTool="./lib/dcpTool.exe",out_path = "./data"):
        
        file_path = os.path.realpath(file_path)
        dcpTool = os.path.realpath(dcpTool)
        out_path = os.path.realpath(out_path)

        print(file_path)
        print(dcpTool)
        print(out_path)

        self.file_path = file_path
        self.dcpTool = dcpTool

        # os.path.basename()获取文件名
        filename = os.path.basename(file_path) 
        # os.path.splitext()获取文件名和拓展名
        filename_without_ext, file_extension = os.path.splitext(filename)

        if file_extension == '.xml':
            self.out_path = file_path
        else:
            self.out_path = os.path.join(out_path,filename_without_ext+".xml")

    def convert2xml(self):
        if self.file_path.endswith('.dcp'):
            # result = subprocess.run([self.dcpTool, '-d' ,self.file_path,self.out_path ], stdout=subprocess.PIPE)
            # print(result.stdout.decode())
            command = f" \"{self.dcpTool}\" -d \"{self.file_path}\" \"{self.out_path}\""
            print("Execute: ")
            print(command)

            subprocess.check_output(command, stderr=subprocess.STDOUT,shell=True) 
            
            if os.path.isfile(self.out_path):
                print("Unpack Success") 
            
            return

        if self.file_path.endswith('.xml'):
            print("Xml File loaded,no need unpack")
        

    def parse(self):
        self.parser = DCPParser(self.out_path)

    
    def save_modified(self,xml_file_outpath,dcp_outpath=None,rename_profileName=None,tags_to_exclude=[]):

        xml_file_outpath = os.path.realpath(xml_file_outpath)
        self.parser.save_modified(xml_file_outpath, rename_profileName= rename_profileName , tags_to_exclude=tags_to_exclude)

         # os.path.basename()获取文件名
        filename = os.path.basename(xml_file_outpath) 
        # os.path.splitext()获取文件名和拓展名
        filename_without_ext, file_extension = os.path.splitext(filename)

        if dcp_outpath is not None:
            dcp_outpath = os.path.realpath(dcp_outpath)
            dcp_file_path = os.path.join(dcp_outpath,filename_without_ext+".dcp") 
            
            command = f" \"{self.dcpTool}\" -c \"{xml_file_outpath}\" \"{dcp_file_path}\""
            print("Execute: ")
            print(command)
            
            subprocess.check_output(command, stderr=subprocess.STDOUT,shell=True) 

            if os.path.isfile(dcp_file_path):
                print("Pack Success") 
            else:
                print("Pack failed")
