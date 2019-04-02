import math
import pandas as pd
from collections import OrderedDict
from pyswmm import Simulation
import pyswmm
import numpy as np
import datetime
import pytz
import influxdb
import random
import math
from affine import Affine

# CLASSES
class swmmINP:
    
    def __init__(self,inpF,*args,**kwargs):
        self.inpF = inpF
        self.args = args
        self.kwargs = kwargs
        self.warnings = []

        # Handle keyword args.
        self.set_kwargs()
            
        self.make_sections()
        # self.prep_dicts() # think this should be handled if you know what headers you have
        
    # SET
    def set_kwargs(self):
        # DEFAULTS
        self.offset = 0.0
        self._min_slope = 0.0
        self.headers = ['[TITLE]','[OPTIONS]','[EVAPORATION]','[RAINGAGES]','[SUBCATCHMENTS]',
           '[SUBAREAS]','[INFILTRATION]','[JUNCTIONS]','[OUTFALLS]','[STORAGE]','[CONDUITS]',
           '[PUMPS]','[ORIFICES]','[WEIRS]','[XSECTIONS]','[LOSSES]','[CONTROLS]','[INFLOWS]',
           '[DWF]','[HYDROGRAPHS]','[RDII]','[CURVES]','[TIMESERIES]','[PATTERNS]','[REPORT]',
           '[TAGS]','[MAP]','[COORDINATES]','[VERTICES]','[Polygons]','[SYMBOLS]','[PROFILES]']
        
        kw = self.kwargs.keys()
        if 'offset' in kw:
            self.offset = self.kwargs['offset']
        
        if 'min_slope' in kw:
            self._min_slope = self.kwargs['slope']

        if 'headers' in kw:
            self.headers = self.kwargs['headers']

    def set_dicts(self):
        self.make_title_dict()
        self.make_xsec_dict()
        self.make_curves_dictionary()
        self.make_conduit_dictionary()
        self.make_junction_dictionary()
        self.make_storage_dictionary()
        self.make_subcatchment_dictionary()
        self.make_subareas_dictionary()
        self.make_outfall_dictionary()
        self.make_orifice_dictionary()
        self.make_pump_dictionary()
        self.make_options_dictionary()
        self.make_raingauges_dict()
        self.make_timeseries_dict()
        self.make_evap_dict()
        self.make_infiltration_dict()
        self.make_losses_dict()
        self.make_controls_dict()
        self.make_inflows_dict()
        self.make_weir_dict()
        self.make_dwf_dict()
        self.make_hydrograph_dict()
        self.make_rdii_dict()
        self.make_report_dict()
        self.make_patterns_dict()
        self.make_tags_dict()
        self.make_profiles_dict()


        self.calc_slope()
        self.calc_qfull()



    # ------ GEOMETRY SECTION ------
    # Includes the affine transformation equations

    def set_geo_dicts(self):
        self.make_map_dict()
        self.make_coords_dict()
        self.make_verts_dict()
        self.make_polygons_dict()
        self.make_symbols_dict()

    def make_anon(self):

        # Perform Affine Transformations on the Geometry and Coordinates of the
        # SWMM input file. Relevant Input sections:
        # [MAP], [COORDINATES], [VERTICES], [POLYGONS], and [SYMBOLS]
        # Uses the Affine module.

        self.rotate('xy','xy_1')

        x_min,x_max,y_min,y_max = self.get_min_max('xy_1')
        
        self.scale(x_min,x_max,y_min,y_max,'xy_1','xy_2')

        x_min,x_max,y_min,y_max = self.get_min_max('xy_2')

        x_trans,y_trans = self.translate(x_min,x_max,y_min,y_max,'xy_2','xy_3')

        self.map['dim_transformed'] = [0,0,1.05 * (x_max + x_trans),1.05 * (y_max + y_trans)]

        self.anon = True

    def rotate(self,start_key,end_key):

        Rotation_Angle = random.random()*360
        A_rot = Affine.rotation(Rotation_Angle)

        # Rotate Nodes:
        for node in self.coords:
            self.coords[node][end_key] = A_rot * self.coords[node][start_key]

         # Rotate Symbols
        for sym in self.symbols:
            self.symbols[sym][end_key] = A_rot * self.symbols[sym][start_key]

        # Rotate Vertices
        for v in self.verts:
            self.verts[v][end_key] = []
            for tup in self.verts[v][start_key]:
                self.verts[v][end_key].append(A_rot * tup)

        # Rotate Polygons
        for p in self.polygons:
            self.polygons[p][end_key] = []
            for tup in self.polygons[p][start_key]:
                self.polygons[p][end_key].append(A_rot * tup)

    def get_min_max(self,key):
        x_min = math.inf
        x_max = -math.inf
        y_min = math.inf
        y_max = -math.inf

        for n in self.coords:
            # Min
            if self.coords[n][key][0] < x_min:
                x_min = self.coords[n][key][0]
            if self.coords[n][key][1] < y_min:
                y_min = self.coords[n][key][1]
                
            # Max
            if self.coords[n][key][0] > x_max:
                x_max = self.coords[n][key][0]
            if self.coords[n][key][1] > y_max:
                y_max = self.coords[n][key][1]

        for s in self.symbols:
            # Min
            if self.symbols[s][key][0] < x_min:
                x_min = self.symbols[s][key][0]
            if self.symbols[s][key][1] < y_min:
                y_min = self.symbols[s][key][1]
                
            # Max
            if self.symbols[s][key][0] > x_max:
                x_max = self.symbols[s][key][0]
            if self.symbols[s][key][1] > y_max:
                y_max = self.symbols[s][key][1]

        for p in self.polygons:
            for tup in self.polygons[p][key]:
                # Min
                if tup[0] < x_min:
                    x_min = tup[0]
                if tup[1] < y_min:
                    y_min = tup[0]
                    
                # Max
                if tup[0] > x_max:
                    x_max = tup[0]
                if tup[1] > y_max:
                    y_max = tup[1]

        for v in self.verts:
            for tup in self.verts[v][key]:
                # Min
                if tup[0] < x_min:
                    x_min = tup[0]
                if tup[1] < y_min:
                    y_min = tup[1]
                    
                # Max
                if tup[0] > x_max:
                    x_max = tup[0]
                if tup[1] > y_max:
                    y_max = tup[1]


        return x_min,x_max,y_min,y_max

    def scale(self,x_min,x_max,y_min,y_max,key_1,key_2):
        x_scale = 100000/(x_max - x_min)
        y_scale = 100000/(y_max - y_min)

        A_scale = Affine.scale(x_scale,y_scale)

        for node in self.coords:
            self.coords[node][key_2] = A_scale * self.coords[node][key_1]
            
        for sym in self.symbols:
            self.symbols[sym][key_2] = A_scale * self.symbols[sym][key_1]

        for p in self.polygons:
            self.polygons[p][key_2] = []
            for tup in self.polygons[p][key_1]:
                self.polygons[p][key_2].append(A_scale * tup)
                
        for v in self.verts:
            self.verts[v][key_2] = []
            for tup in self.verts[v][key_1]:
                self.verts[v][key_2].append(A_scale * tup)

    def translate(self,x_min,x_max,y_min,y_max,key_1,key_2):
        x_trans = -x_min + 0.05 * (x_max - x_min)
        y_trans = -y_min + 0.05 * (y_max - y_min)

        A_trans = Affine.translation(x_trans,y_trans)

        for node in self.coords:
            self.coords[node][key_2] = A_trans * self.coords[node][key_1]
            
        for sym in self.symbols:
            self.symbols[sym][key_2] = A_trans * self.symbols[sym][key_1]

        for p in self.polygons:
            self.polygons[p][key_2] = []
            for tup in self.polygons[p][key_1]:
                self.polygons[p][key_2].append(A_trans * tup)
                
        for v in self.verts:
            self.verts[v][key_2] = []
            for tup in self.verts[v][key_1]:
                self.verts[v][key_2].append(A_trans * tup)

        return x_trans,y_trans



    # ------ MAKE SECTION ------
    # Organize input file into dictionary dataframes
    
    def make_conduit_dictionary(self):
        self.conduits = {}
        for l in self._sections['[CONDUITS]']:
            a = l.split()

            self.conduits[a[0]] = {
                'from_node': a[1],
                'to_node': a[2],
                'length': float(a[3]),
                'roughness': float(a[4]),
                'in_offset': float(a[5]),
                'out_offset': float(a[6]),
                'init_flow': float(a[7]),
                'max_flow': float(a[8]),
            }

        for c in self.conduits:
            self.conduits[c].update(self.xsections[c])

        self.calc_conduit_vol()

    def make_controls_dict(self):
        control_list = [l for l in self._sections['[CONTROLS]']]
        self.controls = {}
        for i in control_list:
            if 'RULE' in i:
                self.controls[i.split()[1]] = []

        for i in control_list:
            if 'RULE' in i:
                rule = i.split()[1]
            else:
                self.controls[rule].append(i)

    def make_coords_dict(self):
        self.coords = OrderedDict()

        for l in self._sections['[COORDINATES]']:
            a = l.split()
            self.coords[a[0]] = {
                'x': float(a[1]),
                'y': float(a[2]),
                'xy': (float(a[1]),float(a[2]))
            }
        
    def make_curves_dictionary(self):
        self.curves = {}

        for l in self._sections['[CURVES]']:
            a = l.split()

            if len(a) == 4:
                self.curves[a[0]] = {
                    'type':a[1],
                    'x_val':[a[2]],
                    'y_val':[a[3]]
                }

            if len(a) == 3:
                self.curves[a[0]]['x_val'].append(a[1])
                self.curves[a[0]]['y_val'].append(a[2])

    def make_dwf_dict(self):
        self.dwf = {}

        for l in self._sections['[DWF]']:
            a = l.split()
            self.dwf[a[0]] = {
                'Constituent' : a[1],
                'Baseline'    : a[2]
            }

            if len(a) > 3:
                self.dwf[a[0]]['Patters'] = a[3:]
                
    def make_evap_dict(self):
        evap_raw = self._sections['[EVAPORATION]']
        self.evaporation = {
            'constant' : float(evap_raw[0].split()[1]),
            'dry_only' : evap_raw[1].split()[1]    
        }
        del evap_raw

    def make_hydrograph_dict(self):
        self.hydrographs = dict()

        for l in self._sections['[HYDROGRAPHS]']:
            a = l.split()
            
            if a[0] not in self.hydrographs.keys():
                self.hydrographs[a[0]] = dict()
                self.hydrographs[a[0]]['Rain Gage'] = a[1]
            
            else:
                self.hydrographs[a[0]]['{0} {1}'.format(a[1],a[2])] = {
                    'R'      : a[3],
                    'T'      : a[4],
                    'K'      : a[5],
                    'Dmax'   : a[6],
                    'Drecov' : a[7],
                    'Dinit'  : a[8]
                }
        
    def make_infiltration_dict(self):
        self.infiltration = {}
        for l in self._sections['[INFILTRATION]']:
            a = l.split()

            self.infiltration[a[0]] = {
                'max_rate' : float(a[1]),
                'min_rate' : float(a[2]),
                'decay' : float(a[3]),
                'dry_time' : float(a[4]),
                'max_infil' : float(a[5])
            }

    def make_inflows_dict(self):
        self.inflows = {}
        for l in self._sections['[INFLOWS]']:
            a = l.split()

            self.inflows[a[0]] = {
                'constituent' : a[1],
                'time_series' : a[2],
                'type' : a[3],
                'm_factor' : float(a[4]),
                's_factor' : float(a[5])
            }

            try:
                a[0]['baseline'] = a[6]
                a[0]['pattern'] = a[7]
            except:
                pass
            
    def make_junction_dictionary(self):
        self.junctions = {}

        for l in self._sections['[JUNCTIONS]']:
            a = l.split()

            self.junctions[a[0]] = {
                'elevation': float(a[1]),
                'max_depth': float(a[2]),
                'init_depth': float(a[3]),
                'sur_depth': float(a[4]),
                'a_ponded': float(a[5])
            }

    def make_losses_dict(self):
        self.losses = {}
        for l in self._sections['[LOSSES]']:
            a = l.split()

            self.losses[a[0]] = {
                'k_entry' : float(a[1]),
                'k_exit' : float(a[2]),
                'k_avg' : float(a[3]),
                'flap_gate' : a[4],
                'seepage' : float(a[5])
            }

    def make_map_dict(self):
        dims = self._sections['[MAP]'][0].split() 
        # I think the dimensions are in the order of East, South, West, North
        self.map = {
            'dim': {
                'west' : float(dims[1]),
                'south' : float(dims[2]),
                'east' : float(dims[3]),
                'north' : float(dims[4])
            },

            'units' : self._sections['[MAP]'][1].split()[1]
        }

    def make_outfall_dictionary(self):
        self.outfalls = {}

        for l in self._sections['[OUTFALLS]']:
            a = l.split()
            self.outfalls[a[0]] = {
                'elevation': float(a[1]),
                'type': a[2],
            }

            if a[2] == 'FREE':
                self.outfalls[a[0]]['stage_data'] = 0.0
                self.outfalls[a[0]]['gated'] = a[3]
            else:
                self.outfalls[a[0]]['stage_data'] = float(a[3]),
                self.outfalls[a[0]]['gated'] = a[4]

    def make_orifice_dictionary(self):
        self.orifices = {}

        for l in self._sections['[ORIFICES]']:
            a = l.split()
            self.orifices[a[0]] = {
                'from_node': str(a[1]),
                'to_node': str(a[2]),
                'type':a[3],
                'offset':float(a[4]),
                'Cd':float(a[5]),
                'gated':a[6],
                'close_time':float(a[7])
            }

        for l in self._sections['[XSECTIONS]']:
            a = l.split()

            try:
                self.orifices[a[0]]['shape'] = a[1]
                self.orifices[a[0]]['geom1'] = float(a[2])
                self.orifices[a[0]]['geom2'] = float(a[3])
                self.orifices[a[0]]['geom3'] = float(a[4])
                self.orifices[a[0]]['geom4'] = float(a[5])
                self.orifices[a[0]]['barrels'] = int(a[6])
                self.orifices[a[0]]['culvert'] = float(a[7])
            except:
                pass

    def make_options_dictionary(self):
        self.options = {}

        for l in self._sections['[OPTIONS]']:
            a = l.split()
            self.options[a[0]]=a[1]

        step_str = self.options['ROUTING_STEP'].split(':')
        timestep_sec = float(step_str[0])*60*60 + float(step_str[1])*60 + float(step_str[2])
        self.options['ROUTING_STEP'] = timestep_sec

    def make_patterns_dict(self):
        self.patterns = dict()

        for l in self._sections['[PATTERNS]']:
            a = l.split()
            
            if a[0] not in self.patterns.keys():
                self.patterns[a[0]] = {
                    'Type'        : a[1],
                    'Multipliers' : [a[2:]]
                }
                
            else:
                self.patterns[a[0]]['Multipliers'].append(a[1:])

    def make_polygons_dict(self):
        self.polygons = OrderedDict()

        for p in self._sections['[Polygons]']:
            l = p.split()
            try:
                self.polygons[l[0]]['x'].append(float(l[1]))
                self.polygons[l[0]]['y'].append(float(l[2]))
                self.polygons[l[0]]['xy'].append((float(l[1]),float(l[2])))
            except:
                self.polygons[l[0]] = {
                    'x' : [float(l[1])],
                    'y' : [float(l[2])],
                    'xy': [(float(l[1]),float(l[2]))]
                }

    def make_profiles_dict(self):
        self.profiles = dict()

    def make_pump_dictionary(self):
        self.pumps = {}
        for l in self._sections['[PUMPS]']:
            a = l.split()
            self.pumps[a[0]] = {
                'from_node':a[1],
                'to_node':a[2],
                'pump_curve':a[3],
                'status':a[4],
                'startup':a[5],
                'shutoff':a[6]
            }

        for p in self.pumps:
            self.pumps[p]['curve_info'] = self.curves[self.pumps[p]['pump_curve']]

    def make_raingauges_dict(self):
        self.raingauges = {}
        for l in self._sections['[RAINGAGES]']:
            a = l.split()
            self.raingauges[a[0]] = {
                'format' : a[1],
                'interval' : float(a[2]),
                'SCF' : float(a[3]),
                'source1' : a[4],
                'source2' : a[5]
            }

    def make_rdii_dict(self):
        self.rdii = dict()

        for l in self._sections['[RDII]']:
            a = l.split()

            self.rdii[a[0]] = {
                'Unit Hydrograph' : a[1],
                'Sewer Area'      : a[2]
            }

    def make_report_dict(self):
        self.report = dict()

        for l in self._sections['[REPORT]']:
            a = l.split()
            self.report[a[0]] = a[1]

    def make_sections(self):
        with open(self.inpF) as f:
            contents = f.read()
        
        self._sections = {}
        for header in self.headers:
            self._sections[header] = contents.find(header)
            
        sort = sorted(self._sections.items(), key=lambda x: x[1])
        
        for i in range(0,len(sort)):
            if i < len(sort)-1:
                a = [sort[i][1],sort[i+1][1]]
            else:
                a = [sort[i][1],len(contents)]
                
            section_content = contents[a[0]:a[1]]
            h = section_content.split('\n')[0]

            self._sections[h] = []
            
            for l in section_content.split('\n'):
                if not l:
                    pass
                elif l[0].isalnum():
                    self._sections[h].append(l)
                else:
                    pass

    def make_storage_dictionary(self):
        self.storages = {}

        for l in self._sections['[STORAGE]']:
            a = l.split()

            self.storages[a[0]] = {
                'elevation':float(a[1]),
                'max_depth':float(a[2]),
                'init_depth':float(a[3]),
                'shape':a[4],
                #Curve Name/Params
                #N/A
                #Fevap
                #PSI
                #Ksat
                #IMD
            }

            if a[4] == 'FUNCTIONAL':
                self.storages[a[0]]['A'] = float(a[5])
                self.storages[a[0]]['B'] = float(a[6])
                self.storages[a[0]]['C'] = float(a[7])

            elif a[4] == 'TABULAR':
                self.storages[a[0]]['curve_name'] = a[5]
                self.storages[a[0]]['curve_info'] = self.curves[a[5]]
                self.storages[a[0]]['curve_info']['x_val'] = [float(i) for i in self.storages[a[0]]['curve_info']['x_val']]
                self.storages[a[0]]['curve_info']['y_val'] = [float(i) for i in self.storages[a[0]]['curve_info']['y_val']]
                self.storages[a[0]]['curve_info']['vol'] = []

            else:
                print(a[0] + ' does not have depth v area info.')

        self.calc_storage_vol()

    def make_subareas_dictionary(self):
        self.subareas = {}

        for l in self._sections['[SUBAREAS]']:
            a = l.split()

            self.subareas[a[0]] = {

                'N-Imperv' : float(a[1]),
                'N-Perv'   : float(a[2]),
                'S-Imperv' : float(a[3]),
                'S-Perv'   : float(a[4]),
                'PctZero'  : float(a[5]),
                'RouteTo'  : a[6],
                'PctRouted': float(a[7])
            }
        
    def make_subcatchment_dictionary(self):
        self.subcatchments = {}

        for l in self._sections['[SUBCATCHMENTS]']:
            a = l.split()

            self.subcatchments[a[0]] = {
                'rain_gage': a[1],
                'outlet': a[2],
                'area': float(a[3]),
                'per_imperv': float(a[4]),
                'width': float(a[5]),
                'slope': float(a[6]),
                'curblen': float(a[7])
            }

    def make_symbols_dict(self):
        self.symbols = OrderedDict()

        for p in self._sections['[SYMBOLS]']:
            l = p.split()

            self.symbols[l[0]] = {
                'x' : float(l[1]),
                'y' : float(l[2]),
                'xy': (float(l[1]),float(l[2]))
            }

    def make_tags_dict(self):
        # this looked empty
        # can't really tell 
        # what to do.
        self.tags = dict()
            
    def make_timeseries_dict(self):
        self.timeseries = {}
        for l in self._sections['[TIMESERIES]']:
            a = l.split()

            if len(a) == 4:
                try:
                    self.timeseries[a[0]]['date'].append(a[1])
                    self.timeseries[a[0]]['time'].append(a[2])
                    self.timeseries[a[0]]['value'].append(a[3])
                except:
                    self.timeseries[a[0]] = {
                        'date' : [a[1]],
                        'time' : [a[2]],
                        'value' : [a[3]]
                    }
            else:
                try:
                    self.timeseries[a[0]]['date'].append('')
                    self.timeseries[a[0]]['time'].append(a[1])
                    self.timeseries[a[0]]['value'].append(a[2])
                except:
                    self.timeseries[a[0]] = {
                        'date' : [''],
                        'time' : [a[1]],
                        'value' : [a[2]]
                    }

    def make_title_dict(self):
        self.title = {}
        self.title['title'] = self._sections['[TITLE]'][0]
        self.title['notes'] = self._sections['[TITLE]'][1]

    def make_verts_dict(self):
        self.verts = OrderedDict()

        for v in self._sections['[VERTICES]']:
            l = v.split()
            try:
                self.verts[l[0]]['x'].append(float(l[1]))
                self.verts[l[0]]['y'].append(float(l[2]))
                self.verts[l[0]]['xy'].append((float(l[1]),float(l[2])))

            except:
                self.verts[l[0]] = {
                    'x' : [float(l[1])],
                    'y' : [float(l[2])],
                    'xy': [(float(l[1]),float(l[2]))]
                }

    def make_weir_dict(self):
        self.weirs = {}

        for l in self._sections['[WEIRS]']:
            a = l.split()

            self.weirs[a[0]] = {
                'From Node' : a[1],
                'To Node'   : a[2],
                'Type'      : a[3],
                'CrestHt'   : float(a[4]),
                'Qcoeff'    : float(a[5]),
                'Gated'     : a[6],
                'EndCon'    : float(a[7]),
                'EndCoeff'  : float(a[8]),
                'Surcharge' : a[9],
            }

            try:
                self.weirs[a[0]]['RoadWidth'] = a[10]
                self.weirs[a[0]]['RoadSurf']  = a[11]
            except:
                pass

    def make_xsec_dict(self):
        self.xsections = {}
        for l in self._sections['[XSECTIONS]']:
            a = l.split()

            if len(a) != 8 and len(a) != 6:
                self.warnings.append(a[0] + ' ' + str(len(a)))

            self.xsections[a[0]] = {
                'shape':a[1],
                'geom1' : float(a[2]),
                'geom2' : float(a[3]),
                'geom3' : float(a[4]),
                'geom4' : float(a[5]),
            }

            try:
                self.xsections[a[0]]['barrels'] = int(a[6])
                self.xsections[a[0]]['culvert'] = float(a[7])
            except:
                self.xsections[a[0]]['barrels'] = 1
                self.xsections[a[0]]['culvert'] = 0

        self.calc_xsec_area()

    

    # CALC
    def calc_xsec_area(self):
        for item in self.xsections:
            if self.xsections[item]['shape'] == 'CIRCULAR':
                self.xsections[item]['area'] = ( self.xsections[item]['geom1'] / 2 ) ** 2 * math.pi * self.xsections[item]['barrels']
            elif self.xsections[item]['shape'] == 'RECT_CLOSED' or self.xsections[item]['shape'] == 'RECT_OPEN':
                self.xsections[item]['area'] = self.xsections[item]['geom1'] * self.xsections[item]['geom2'] * self.xsections[item]['barrels']
            elif self.xsections[item]['shape'] == 'TRIANGULAR':
                self.xsections[item]['area'] = 0.5 * self.xsections[item]['geom1'] * self.xsections[item]['geom2'] * self.xsections[item]['barrels']
            else:
                self.warnings.append(item + ' ' + self.xsections[item]['shape'] + ' not yet calculated area')
                self.xsections[item]['area'] = 1.0

    def convert(self, conversion):
        to_convert = [self.coords,self.verts,self.polygons]
        for i in to_convert:
            for k in i:
                x = np.array(i[k]['x']) / conversion
                y = np.array(i[k]['y']) / conversion

                i[k]['x'] = x.tolist()
                i[k]['y'] = y.tolist()

    def calc_datum_conversion(self,key_name):
        # apply offset to variables that have elevation:
        #     - Junctions
        #     - Storages
        #     - Outfalls
        for point in self.junctions:
            self.junctions[point][key_name] = self.junctions[point]['elevation'] - self.offset
        for point in self.storages:
            self.storages[point][key_name] = self.storages[point]['elevation'] - self.offset
        for point in self.outfalls:
            self.outfalls[point][key_name] = self.outfalls[point]['elevation'] - self.offset            
            
    def calc_slope(self):
        for item in self.conduits:
            if self.conduits[item]['from_node'] in self.junctions.keys():
                e1 = self.junctions[self.conduits[item]['from_node']]['elevation']+self.conduits[item]['in_offset']
            elif self.conduits[item]['from_node'] in self.storages.keys():
                e1 = self.storages[self.conduits[item]['from_node']]['elevation']+self.conduits[item]['in_offset']
            else:
                e1 = 1

            if self.conduits[item]['to_node'] in self.junctions.keys():
                e2 = self.junctions[self.conduits[item]['to_node']]['elevation']+self.conduits[item]['out_offset']
            elif self.conduits[item]['to_node'] in self.storages.keys():
                e2 = self.storages[self.conduits[item]['to_node']]['elevation']+self.conduits[item]['out_offset']
            else:
                e2 = 1

            if e1==1 or e2==1:
                self.conduits[item]['slope_flag'] = True
            else:
                self.conduits[item]['slope_flag'] = False

            slope = (e1 - e2)/self.conduits[item]['length']

            if slope < self._min_slope:
                slope = self._min_slope
                self.conduits[item]['slope_flag'] = True

            self.conduits[item]['slope'] = slope
            
    def calc_qfull(self):
        for item in self.conduits:
            if self.conduits[item]['shape'] == 'CIRCULAR':
                # compute Qfull as pipe full manning equation
                self.conduits[item]['q_full'] = (self.conduits[item]['geom1']**(8/3)*self.conduits[item]['slope']**(1/2))/(4**(5/3)*self.conduits[item]['roughness'])*math.pi
            elif self.conduits[item]['shape'] == 'RECT_CLOSED':
                # Compute q_full as manning equation of pipe with manning eq with depth as 0.95
                self.conduits[item]['q_full'] = (1.49/self.conduits[item]['roughness']) * (0.95 * self.conduits[item]['geom1'] * self.conduits[item]['geom2']) * (self.conduits[item]['geom2'] * 0.95 * self.conduits[item]['geom1'] / (self.conduits[item]['geom2'] + 2 * 0.95 * self.conduits[item]['geom1']))**(2/3)
            else:
                self.conduits[item]['q_full'] = 1;

    def calc_conduit_vol(self):
        for element in self.conduits:
            self.conduits[element]['vol'] = self.conduits[element]['area'] * self.conduits[element]['length']

    def calc_storage_vol(self):
        for element in self.storages:
            if self.storages[element]['shape'] == 'FUNCTIONAL':
                self.storages[element]['total_storage'] = ( self.storages[element]['A'] * self.storages[element]['max_depth']**(self.storages[element]['B'] + 1) / ( self.storages[element]['B'] + 1 ) ) + ( self.storages[element]['C'] * self.storages[element]['max_depth'] )
            elif self.storages[element]['shape'] == 'TABULAR':
                '''
                'a' is any generic area under a curve.
                FROM SWMM SOURCE CODE:
                    The area within each interval i of the table is given by:
                    Integral{ y(x)*dx } from x(i) to x
                    where y(x) = y(i) + s*dx
                    dx = x - x(i)
                    s = [y(i+1) - y(i)] / [x(i+1) - x(i)]
                    This results in the following expression for a(i):
                    a(i) = y(i)*dx + s*dx*dx/2
                '''
                x = self.storages[element]['curve_info']['x_val']
                y = self.storages[element]['curve_info']['y_val']
                v = [0.0]
                for i in range(1,len(x)):
                    h = x[i] - x[i-1]
                    a = (y[i-1] + y[i])*h/2

                    v.append(v[-1]+a)

                self.storages[element]['curve_info']['vol'] = v

                '''Generally speaking, if the storage curve is tabular then the last value seems to be a really large number that would more closely align with flooding of a "storage" element that is really maybe just a manhole. So, we'll call the total storage approximately being equal to the second to last volume value.'''
                self.storages[element]['total_storage'] = self.storages[element]['curve_info']['vol'][-2]
            else:
                pass


    # ------ MODEL OUTPUT ------

    def save_input(self, **kwargs):

        filename = self.inpF

        if 'filename' in kwargs.keys():
            filename = kwargs['filename']

        self.file_obj = open(filename,'w')


        # Call Write sections Here
        self.write_title()
        self.write_options()
        self.write_evaporation()
        self.write_raingages()
        self.write_subcatchments()
        self.write_subareas()
        self.write_infiltration()
        self.write_junctions()
        self.write_outfalls()
        self.write_storages()
        self.write_conduits()
        self.write_pumps()
        self.write_orifices()
        self.write_weirs()
        self.write_xsections()
        self.write_losses()
        self.write_controls()
        self.write_inflows()
        self.write_dwf()
        self.write_hydrographs()
        self.write_rdii()
        self.write_curves()
        self.write_timeseries()
        self.write_patterns()
        self.write_report()
        self.write_tags()
        self.write_map()
        self.write_coordinates()
        self.write_vertices()
        self.write_polygons()
        self.write_symbols()
        self.write_profiles()

        self.file_obj.close()


    def write_title(self):
        self.file_obj.write('[TITLE]\n')
        self.file_obj.write(';;Project Title/Notes\n')
        self.file_obj.write(self.title['title'] + '\n')
        self.file_obj.write(self.title['notes'] + '\n')
        self.file_obj.write('\n')

    def write_options(self):
        self.file_obj.write('[OPTIONS]\n')
        self.file_obj.write(';;Option \t \t Value\n')

        for k in self.options:
            self.file_obj.write('{0}\t{1}\n'.format(k,self.options[k]))

            if k == 'ROUTING_STEP' or k == 'SKIP_STEADY_STATE' or k == 'THREADS':
                self.file_obj.write('\n')

    def write_evaporation(self):
        self.file_obj.write('[EVAPORATION]\n')
        self.file_obj.write(';;Data Source    Parametes\n')
        self.file_obj.write(';;-------------- -----------------\n')

        for e in self.evaporation:
            self.file_obj.write('{0}\t{1}\n'.format(e,self.evaporation[e]))

        self.file_obj.write('\n')

    def write_raingages(self):
        self.file_obj.write('[RAINGAGES]\n')
        self.file_obj.write(';;Name           Format    Interval SCF      Source \n')
        self.file_obj.write(';;-------------- --------- ------ ------ ----------\n')

        for r in self.raingauges:
            line = '{0}\t\t\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                r,
                str(self.raingauges[r]['format']),
                str(self.raingauges[r]['interval']),
                str(self.raingauges[r]['SCF']),
                str(self.raingauges[r]['source1']),
                str(self.raingauges[r]['source2']),
                )

            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_subcatchments(self):
        self.file_obj.write('[SUBCATCHMENTS]\n')
        self.file_obj.write(';;Name           Rain Gage        Outlet           Area     %Imperv  Width    %Slope   CurbLen  SnowPack\n')
        self.file_obj.write(';;-------------- ---------------- ---------------- -------- -------- -------- -------- -------- ----------------\n')
        
        for sub in self.subcatchments:
            s = self.subcatchments[sub]
            self.file_obj.write(sub + '\t\t\t' + s['rain_gage'] + '\t\t' + s['outlet'] + '\t\t' + str(s['area']) + '\t' + str(s['per_imperv']) + '\t' + str(s['width']) + '\t' + str(s['slope']) + '\t' + str(s['curblen']) + '\n')
        self.file_obj.write('\n')

    def write_subareas(self):
        self.file_obj.write('[SUBAREAS]\n')
        self.file_obj.write(';;Subcatchment   N-Imperv   N-Perv     S-Imperv   S-Perv     PctZero    RouteTo    PctRouted \n')
        self.file_obj.write(';;-------------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n')

        for sub in self.subareas:
            s = self.subareas[sub]
            line = '{0} \t\t {1}   {2}   {3}   {4}   {5}  {6}   {7}  \n'.format(
                sub,
                str(s['N-Imperv']),
                str(s['N-Perv']),
                str(s['S-Imperv']),
                str(s['S-Perv']),
                str(s['PctZero']),
                str(s['RouteTo']),
                str(s['PctRouted'])
                )

            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_infiltration(self):
        self.file_obj.write('[INFILTRATION]\n')
        self.file_obj.write(';;Subcatchment   MaxRate    MinRate    Decay      DryTime    MaxInfil  \n')
        self.file_obj.write(';;-------------- ---------- ---------- ---------- ---------- ----------\n')

        for sub in self.infiltration:
            i = self.infiltration[sub]

            line = '{0} \t\t\t {1}   {2}   {3}   {4}   {5} \n'.format(
                sub,
                str(i['max_rate']),
                str(i['min_rate']),
                str(i['decay']),
                str(i['dry_time']),
                str(i['max_infil'])
                )

            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_junctions(self):
        self.file_obj.write('[JUNCTIONS]\n')
        self.file_obj.write(';;Name           Elevation  MaxDepth   InitDepth  SurDepth   Aponded   \n')
        self.file_obj.write(';;-------------- ---------- ---------- ---------- ---------- ----------\n')

        for node in self.junctions:
            n = self.junctions[node]

            line = '{0} \t\t\t {1} {2}  {3}   {4}   {5} \n'.format(
                node,
                str(n['elevation']),
                str(n['max_depth']),
                str(n['init_depth']),
                str(n['sur_depth']),
                str(n['a_ponded'])
                )

            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_outfalls(self):
        self.file_obj.write('[OUTFALLS]\n')
        self.file_obj.write(';;Name           Elevation  Type       Stage Data       Gated    Route To        \n')
        self.file_obj.write(';;-------------- ---------- ---------- ---------------- -------- ----------------\n')

        for node in self.outfalls:
            n = self.outfalls[node]

            line = '{0} \t\t\t {1}\t\t {2}'.format(
                node,
                str(n['elevation']),
                n['type']
                )

            if n['type'] == 'FIXED':
                line = line + '\t\t{0}\t\t\t{1}'.format(
                    str(n['stage_data'][0]),
                    n['gated']
                    )
            else:
                line = line + '\t\t\t\t\t\t{0}'.format(
                    n['gated']
                    )


            if 'route_to' in n.keys():
                # do something in here.
                # don't know yet.
                pass


            self.file_obj.write(line+'\n')

        self.file_obj.write('\n')

    def write_storages(self):
        # NOTE: My sample inp file that I am basing this from doesn't include evaporation or infiltration
        # through the storage unit.

        self.file_obj.write('[STORAGE\n')
        self.file_obj.write(';;Name           Elev.    MaxDepth   InitDepth  Shape      Curve Name/Params            N/A      Fevap    Psi      Ksat     IMD     \n')
        self.file_obj.write(';;-------------- -------- ---------- ----------- ---------- ---------------------------- -------- --------          -------- --------\n')

        for stor in self.storages:
            s = self.storages[stor]

            line = '{0} \t\t\t {1} {2}  {3}   {4}'.format(
                stor,
                str(s['elevation']),
                str(s['max_depth']),
                str(s['init_depth']),
                s['shape']
                )

            if s['shape'] == 'FUNCTIONAL':
                line = line + ' {0}  {1}  {2}  0.000000 0.000000\n'.format(
                    str(s['A']),
                    str(s['B']),
                    str(s['C'])
                    )
            else:
                # Is tabular. Insert Curve name.
                line = line +'\t{0}                       0.000000 0.000000\n'.format(
                    s['curve_name']
                    )

            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_conduits(self):
        self.file_obj.write('[CONDUITS]\n')
        self.file_obj.write(';;Name           From Node        To Node          Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow   \n')
        self.file_obj.write(';;-------------- ---------------- ---------------- ---------- ---------- ---------- ---------- ---------- ----------\n')

        for con in self.conduits:
            c = self.conduits[con]
            line = '{0} \t\t\t {1} \t\t\t {2} \t\t\t {3} {4} {5} {6} {7} {8}\n'.format(
                con,
                c['from_node'],
                c['to_node'],
                str(c['length']),
                str(c['roughness']),
                str(c['in_offset']),
                str(c['out_offset']),
                str(c['init_flow']),
                str(c['max_flow'])
                )

            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_pumps(self):
        self.file_obj.write('[PUMPS]\n')
        self.file_obj.write(';;Name           From Node        To Node          Pump Curve       Status   Sartup Shutoff \n')
        self.file_obj.write(';;-------------- ---------------- ---------------- ---------------- ------ -------- --------\n')

        for pump in self.pumps:
            p = self.pumps[pump]

            line = '{0} \t\t\t {1} \t\t {2} \t\t {3} \t {4} \t {5}  {6}\n'.format(
                pump,
                p['from_node'],
                p['to_node'],
                p['pump_curve'],
                p['status'],
                p['startup'],
                p['shutoff']
                )
            self.file_obj.write(line)

        self.file_obj.write('\n')

    def write_orifices(self):
        self.file_obj.write('[ORIFICES]\n')
        self.file_obj.write(';;Name           From Node        To Node          Type         Offset     Qcoeff     Gated    CloseTime \n')
        self.file_obj.write(';;-------------- ---------------- ---------------- ------------ ---------- ---------- -------- ----------\n')

        for orif in self.orifices:
            o = self.orifices[orif]

            line = '{0} \t\t\t {1} \t\t {2} \t\t {3} \t\t {4} \t\t {5} \t\t {6} \t\t {7}\n'.format(
                orif,
                o['from_node'],
                o['to_node'],
                o['type'],
                str(o['offset']),
                str(o['Cd']),
                o['gated'],
                str(o['close_time'])
                )

            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_weirs(self):
        # NOTE: sample input file only had "transverse" type weir
        # Also, no 'RoadWidth' or 'RoadSurf'

        self.file_obj.write('[WEIRS]\n')
        self.file_obj.write(';;Name           From Node        To Node          Type         CrestHt    Qcoeff     Gated    EndCon   EndCoeff   Surcharge  RoadWidth  RoadSurf  \n')
        self.file_obj.write(';;-------------- ---------------- ---------------- ------------ ---------- ---------- -------- -------- ---------- ---------- ---------- ----------\n')

        for we in self.weirs:
            w = self.weirs[we]

            line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(
                we,
                w['From Node'],
                w['To Node'],
                w['Type'],
                str(w['CrestHt']),
                str(w['Qcoeff']),
                w['Gated'],
                str(w['EndCon']),
                str(w['EndCoeff']),
                w['Surcharge']
                )

            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_xsections(self):
        self.file_obj.write('[XSECTIONS]\n')
        self.file_obj.write(';;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels    Culvert   \n')
        self.file_obj.write(';;-------------- ------------ ---------------- ---------- ---------- ---------- ---------- ----------\n')

        for xs in self.xsections:
            x = self.xsections[xs]

            line = '{0}\t\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{7}\n'.format(
                xs,
                x['shape'],
                str(x['geom1']),
                str(x['geom2']),
                str(x['geom3']),
                str(x['geom4']),
                str(x['barrels']),
                str(x['culvert'])
                )

            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_losses(self):
        self.file_obj.write('[LOSSES]\n')
        self.file_obj.write(';;Link           Kentry     Kexit      Kavg       Flap Gate  Seepage   \n')
        self.file_obj.write(';;-------------- ---------- ---------- ---------- ---------- ----------\n')

        for loss in self.losses:
            l = self.losses[loss]

            line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                loss,
                str(l['k_entry']),
                str(l['k_exit']),
                str(l['k_avg']),
                l['flap_gate'],
                str(l['seepage'])
                )
            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_controls(self):
        self.file_obj.write('[CONTROLS]\n')

        for c in self.controls:
            self.file_obj.write('RULE {0}\n'.format(c))
            for l in self.controls[c]:
                self.file_obj.write(l+'\n')
        self.file_obj.write('\n')

    def write_inflows(self):
        # NOTE: NEED TO ADD BASELINE AND PATTERN

        self.file_obj.write('[INFLOWS]\n')
        self.file_obj.write(';;Node           Constituent      Time Series      Type     Mfactor  Sfactor  Baseline Pattern\n')
        self.file_obj.write(';;-------------- ---------------- ---------------- -------- -------- -------- -------- --------\n')

        for inf in self.inflows:
            i = self.inflows[inf]

            line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                inf,
                i['constituent'],
                i['time_series'],
                i['type'],
                str(i['m_factor']),
                str(i['s_factor'])
                )
            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_dwf(self):
        self.file_obj.write('[DWF]\n')
        self.file_obj.write(';;Node           Constituent      Baseline   Patterns  \n')
        self.file_obj.write(';;-------------- ---------------- ---------- ----------\n')

        for dwf in self.dwf:
            d = self.dwf[dwf]
            line = '{0}\t\t{1}\t\t{2}\n'.format(
                dwf,
                d['Constituent'],
                str(d['Baseline'])
                )
            self.file_obj.write(line)
        self.file_obj.write('\n')
        
    def write_hydrographs(self):
        self.file_obj.write('[HYDROGRAPHS]\n')
        self.file_obj.write(';;Hydrograph     Rain Gage/Month  Response R        T        K        Dmax     Drecov   Dinit   \n')
        self.file_obj.write(';;-------------- ---------------- -------- -------- -------- -------- -------- -------- --------\n')

        for name in self.hydrographs:
            self.file_obj.write('{0} {1}\n'.format(name, self.hydrographs[name]['Rain Gage']))
            for key in self.hydrographs[name]:
                if key == 'Rain Gage':
                    pass
                else:
                    h = self.hydrographs[name][key]
                    Month,Response = key.split()
                    line = '{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(
                        name,
                        Month,
                        Response,
                        h['R'],
                        h['T'],
                        h['K'],
                        h['Dmax'],
                        h['Drecov'],
                        h['Dinit']
                    )
                    self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_rdii(self):
        self.file_obj.write('[RDII]\n')
        self.file_obj.write(';;Node           Unit Hydrograph  Sewer Area\n')
        self.file_obj.write(';;-------------- ---------------- ----------\n')

        for node in self.rdii:
            n = self.rdii[node]

            line = '{0}\t{1}\t{2}'.format(
                node,
                n['Unit Hydrograph'],
                n['Sewer Area']
                )
            self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_curves(self):
        self.file_obj.write('[CURVES]\n')
        self.file_obj.write(';;Name           Type       X-Value    Y-Value   \n')
        self.file_obj.write(';;-------------- ---------- ---------- ----------\n')

        for curve in self.curves:
            c = self.curves[curve]

            first_line = '{0} {1} \t {2} {3}\n'.format(
                curve,
                c['type'],
                c['x_val'][0],
                c['y_val'][0]
                )
            self.file_obj.write(first_line)

            for i in range(1,len(c['x_val'])):
                line = '{0} \t\t {1} {2}\n'.format(
                    curve,
                    c['x_val'][i],
                    c['y_val'][i]
                    )
                self.file_obj.write(line)
        self.file_obj.write('\n')

    def write_timeseries(self):
        self.file_obj.write('[TIMESERIES]\n')
        self.file_obj.write(';;Name           Date       Time       Value     \n')
        self.file_obj.write(';;-------------- ---------- ---------- ----------\n')

        for ts in self.timeseries:
            wrt = self.timeseries[ts]
            for d,t,v in zip(wrt['date'],wrt['time'],wrt['value']):
                self.file_obj.write(ts + '\t \t' + str(d) + '\t \t' + str(t) + '\t \t' + str(v) + '\n')
            self.file_obj.write(';\n')
        self.file_obj.write('\n')

    def write_patterns(self):
        self.file_obj.write('[PATTERNS]\n')
        self.file_obj.write(';;Name           Type       Multipliers\n')
        self.file_obj.write(';;-------------- ---------- -----------\n')

        for name in self.patterns:
            p = self.patterns[name]

            line = '{0} \t{1} \t'.format(name,p['Type'])
            for i in p['Multipliers'][0]:
                line += '{0} '.format(i) 
            self.file_obj.write(line + '\n')

            for i in range(1,len(p['Multipliers'])):
                line = name + '\t\t\t'
                for j in p['Multipliers'][i]:
                    line += j + ' '

                self.file_obj.write(line + '\n')
        self.file_obj.write('\n')

    def write_report(self):
        self.file_obj.write('[REPORT]\n')
        self.file_obj.write(';;Reporting Options\n')

        for r in self.report:
            self.file_obj.write(r + ' ' + self.report[r] + '\n')
        self.file_obj.write('\n')

    def write_tags(self):
        self.file_obj.write('[TAGS]\n')
        if self.tags.keys():
            print("Warning: No Tags Added. Contribute to write_tags() in swmmINP class.")
        else:
            pass
        self.file_obj.write('\n')

    def write_map(self):
        self.file_obj.write('[MAP]\n')
        
        line = 'DIMENSIONS '
        # if geometry has been anonymized,
        # use anonymous dimensions.
        if self.anon:
            for i in self.map['dim_transformed']:
                line += str(i) + ' '
        else:
            line += '{0} {1} {2} {3}'.format(
                self.map['dim']['west'],
                self.map['dim']['south'],
                self.map['dim']['east'],
                self.map['dim']['north']
                )
        self.file_obj.write(line + '\n')

        # Write units
        line = 'Units ' + self.map['units'] + '\n'
        self.file_obj.write(line + '\n')
        self.file_obj.write('\n')

    def write_coordinates(self):
        self.file_obj.write('[COORDINATES]\n')
        self.file_obj.write(';;Node           X-Coord            Y-Coord           \n')
        self.file_obj.write(';;-------------- ------------------ ------------------\n')

        if self.anon:
            xy = 'xy_3'
        else:
            xy = 'xy'

        for n in self.coords:
            wrt = '{0}\t\t{1}\t\t{2}\n'.format(
                n,
                self.coords[n][xy][0],
                self.coords[n][xy][1]
            )
            self.file_obj.write(wrt)
        self.file_obj.write('\n')

    def write_vertices(self):
        if self.anon:
            xy = 'xy_3'
        else:
            xy = 'xy'

        self.file_obj.write('[VERTICES]\n')
        self.file_obj.write(';;Link           X-Coord            Y-Coord           \n')
        self.file_obj.write(';;-------------- ------------------ ------------------\n')

        for v in self.verts:
            for xys in self.verts[v][xy]:
                wrt = '{0}\t\t{1}\t\t{2}\n'.format(
                    v,
                    xys[0],
                    xys[1]
                )
                self.file_obj.write(wrt)
        self.file_obj.write('\n')

    def write_polygons(self):
        if self.anon:
            xy = 'xy_3'
        else:
            xy = 'xy'

        self.file_obj.write('[POLYGONS]\n')
        self.file_obj.write(';;Subcatchment   X-Coord            Y-Coord           \n')
        self.file_obj.write(';;-------------- ------------------ ------------------\n')

        for p in self.polygons:
            for xys in self.polygons[p][xy]:
                wrt = '{0}\t\t{1}\t\t{2}\n'.format(
                    p,
                    xys[0],
                    xys[1]
                )
                self.file_obj.write(wrt)
        self.file_obj.write('\n')

    def write_symbols(self):
        if self.anon:
            xy = 'xy_3'
        else:
            xy = 'xy'

        self.file_obj.write('[SYMBOLS]\n')
        self.file_obj.write(';;Gage           X-Coord            Y-Coord           \n')
        self.file_obj.write(';;-------------- ------------------ ------------------\n')
        
        for s in self.symbols:
            wrt = '{0}\t\t{1}\t\t{2}\n'.format(
                s,
                self.symbols[s][xy][0],
                self.symbols[s][xy][1]
            )
            self.file_obj.write(wrt)
        self.file_obj.write('\n')

    def write_profiles(self):
        self.file_obj.write('[PROFILES]\n')
        if self.tags.keys():
            print("Warning: No Profiles Added. Contribute to write_profiles() in swmmINP class.")
        else:
            pass
        self.file_obj.write('\n')


class system():

    def __init__(self, sim, *args, **kwargs):
        self.units = sim.system_units
        self.flood_count = 0.0
        self.tot_flow = 0.0
        self.offset = 0.0
        self.control = True
        self.kwargs = kwargs
        self.group_lookup = dict()
        self.actions = []
        self.control_step = 0.0

        if self.units == 'US':
            self.g = 32.2 # gravity
        else:
            self.g = 9.81 # gravity

        print(self.kwargs)
        kw = self.kwargs.keys()
        if 'offset' in kw:
            self.offset = self.kwargs['offset']
        if 'control' in kw:
            self.control = self.kwargs['control']
        if 'control_step' in kw:
            self.control_step = self.kwargs['control_step']


class ControlPoint:
    def __init__(self,line):
        self.c_name = line[0]
        self.c_type = line[1]
        self.action = float(line[2])
        self.u_name = line[3]
        self.u_type = line[4]
        self.u_param = float(line[5])   # Volume parameter weight
        self.ds_param = float(line[10]) # derivative parameter weight
        self.measure = line[6]
        
        self.location = line[8]
        self.group = int(line[9])
        
        self.recommendations = []
        self.flood_el = float(line[7])
        self.flooding = False
        self.flood_count = 0.0

        self.q_goal = 0.0

    def __str__(self):
        return self.location


    def set_vars(self,nodes,links):
        # Give nodes and links, sets attributes
        # control element
        self.c_var = links[self.c_name]
    
        # upstream element
        if self.u_type == 'junction' or self.u_type == 'storage':
            self.u_var = nodes[self.u_name]
        elif self.u_type == 'link':
            self.u_var = links[self.u_name]

    def get_target_setting(self,run,nodes,links):
        # Could play with the idea of making nodes a global variable.
        control_connects = self.c_var.connections
        upstream = nodes[control_connects[0]]
        downstream = nodes[control_connects[1]]

        h1 = upstream.depth + upstream.invert_elevation
        h2 = downstream.depth + downstream.invert_elevation
        
        current_setting = self.c_var.current_setting # current_setting == hcrown

        pump = False   
        if self.c_type == 'pump':
            pump = True

        if not pump:
            current_height = current_setting * self.cmi['geom1'] # current_height == hcrown
            h_midpt = (current_height / 2) + (upstream.invert_elevation + downstream.invert_elevation) / 2
            hcrest = upstream.invert_elevation + self.cmi['offset']
            
            # inlet submergence
            if h1 < current_height:
                f = (h1 - hcrest) / (current_height - hcrest) # weir equation
            else:
                f = 1.0 # submerged.
    
            # which head to use
            if f < 1.0:
                H = h1 - hcrest
            elif h2 < h_midpt:
                H = h1 - h_midpt
            else:
                H = h1 - h2
            
            # USE CALCULATED HEAD AND DESIRED FLOW TO DETERMINE GATE OPENING ACTION
            
            # no head at orifice
            if H < 0.1 or f <= 0.0:
                self.action = 0.0
                # print('Head too small')
            elif h2 > h1:
                self.action = 0.0
                print('Backward flow condition, orifice closed')
            
            # Weir Flow
            elif (f < 1.0 and H > 0.1):
                A_open = self.q_goal / ( self.cmi['Cd'] * np.sqrt(2*run.g*H) * (2.0/3.0) )
                
                if self.cmi['shape'] == 'CIRCULAR':
                    print(self.c_name)
                    print("Circular does not work yet. Action = 0.0")
                    self.action = 0.0
                else:
                    A_ratio = A_open / ( self.cmi['geom1'] * self.cmi['geom2'] )
                    self.action = A_ratio
            
            # True orifice flow
            else:
                # since q = Cd * A_open * sqrt( 2 g H )
                A_open = self.q_goal / ( self.cmi['Cd'] * np.sqrt(2*run.g*H) )
                
                if self.cmi['shape'] == 'CIRCULAR':
                    print(self.c_name)
                    print("Circular does not work yet. Action = 0.0")
                    self.action = 0.0
                else:
                    A_ratio = A_open / ( self.cmi['geom1'] * self.cmi['geom2'] )
                    self.action = A_ratio
                    

        # Pump is true
        else: 
            if self.cmi['curve_info']['type'] == 'PUMP1':
                print(self.c_name, 'Pump type 1...')

            elif self.cmi['curve_info']['type'] == 'PUMP2':
                # q_out is a function of depth in wet well.
                
                # get q_full from depth
                depth = upstream.depth
                
                q_full = 0.0
                
                n = len(self.cmi['curve_info']['x_val'])
                for  l in range(n-1,-1,-1):
                    # first index in the list. Can't find negative index.
                    if l  == 0:
                        if depth < self.cmi['curve_info']['x_val'][l]:
                            q_full = self.cmi['curve_info']['y_val'][l]

                    else:
                        if depth < self.cmi['curve_info']['x_val'][l] and depth > self.cmi['curve_info']['x_val'][l-1]:
                            q_full = self.cmi['curve_info']['y_val'][l]        

            elif self.cmi['curve_info']['type'] == 'PUMP3':
                head = h2 - h1 # pump pushes water from low head (h1) to higher head (h2)
                
                # x: Head
                # y: Flow
                # x = [float(j) for j in self.cmi['curve_info']['x_val']] 
                # y = [float(j) for j in self.cmi['curve_info']['y_val']]
                
                if head > max(self.cmi['curve_info']['x_val']):
                    head = max(self.cmi['curve_info']['x_val'])
                elif head < min(self.cmi['curve_info']['x_val']):
                    head = min(self.cmi['curve_info']['x_val'])

                # calculate q_full at given head.
                q_full = np.interp(
                    head,
                    np.array(self.cmi['curve_info']['x_val']),
                    np.array(self.cmi['curve_info']['y_val']),
                    )

            elif self.cmi['curve_info']['type'] == 'PUMP4':
                print(self.c_name, 'Pump type 4...')
            
            if self.q_goal == 0.0:
                self.action = 0.0
            else:
                self.action = self.q_goal / q_full 
        
        # if target setting greater than 1, only open to 1.
        self.action = min(max(self.action,0.0),1.0)

        self.check_flooding(run,nodes,links)
        # print(self.c_name,self.c_var.target_setting)
        self.c_var.target_setting = self.action

    def check_flooding(self,run,nodes,links):
        self.flooding = False

        if isinstance(self.u_var, pyswmm.nodes.Node):
            # get node elevation
            elev = self.u_var.depth + self.u_var.invert_elevation
        elif isinstance(self.u_var, pyswmm.links.Link):
            # get inlet node's elevation for the link
            elev =  nodes[self.u_var.inlet_node].depth + nodes[self.u_var.inlet_node].invert_elevation
        
        if elev + run.offset > self.flood_el:
            self.flooding = True
            self.action = 1.0
            self.flood_count = self.flood_count + 1.0

    def get_model_info(self,model):
        # print(self.c_name,self.c_type)
        # print(self.u_name,self.u_type)

        if self.c_type == 'pump':
            self.cmi = model.pumps[self.c_name]
            self.cmi['curve_info']['x_val'] = [float(i) for i in self.cmi['curve_info']['x_val']]
            self.cmi['curve_info']['y_val'] = [float(i) for i in self.cmi['curve_info']['y_val']]
        elif self.c_type == 'orifice':
            self.cmi = model.orifices[self.c_name]

        if self.u_type == 'storage':
            self.umi = model.storages[self.u_name]
            self.max_depth = self.umi['max_depth']
            self.max_vol = self.umi['total_storage']
        elif self.u_type == 'link':
            self.umi = model.conduits[self.u_name]
            self.max_depth = self.umi['geom1']
            self.max_vol = self.umi['vol']

    def get_measure(self):
        if self.measure == 'depth':
            self.now = self.u_var.depth
        elif self.measure == 'flow':
            print(self.name, 'tryna get flow...')
        else:
            print(self.name, 'did not grab measure. Inspect')
    
    def rec_list(self,time_str):
        rec_str = 'REC,site='+self.c_name+' '+'value='+str(self.action) + ' ' + time_str
        self.recommendations.append(rec_str)
        

class DownstreamPoint:
    def __init__(self,line):
        self.d_name = line[0]
        self.d_type = line[1]
        self.measure = line[2]
        self.epsilon = float(line[3])
        self.gamma = float(line[4])
        self.max_depth = float(line[5])
        self.set_point = float(line[6])
        self.set_derivative = float(line[7])
        self.location = line[8]
        self.group = int(line[9])
        self.derivative = 0.0

    def set_vars(self,nodes,links):
        print(self.d_type)
        if self.d_type == 'storage':
            self.d_var = nodes[self.d_name]
        elif self.d_type == 'link':
            self.d_var = links[self.d_name]

    def get_model_info(self,model):
        print('Downstream Points')
        print(self.d_name,self.d_type)

        if self.d_type == 'storage':
            self.dmi = model.storages[self.d_name]
            self.max_vol = self.dmi['total_storage']
        elif self.d_type == 'link':
            self.dmi = model.conduits[self.d_name]
            self.max_vol = self.dmi['vol']

    def get_measure(self):
        if self.measure == 'depth':
            self.now = self.d_var.depth
        elif self.measure == 'flow':
            print('did not write flow yet. Turn around...')
        else:
            print(self.d_name,' did not grab measure. Inspect')


# MAKE
def make_control_points(fn):
    ControlPoints = []
    with open(fn,'r') as f:
        next(f)
        for line in f:
            line = line.strip('\n').split(',')
            ControlPoints.append(ControlPoint(line))
            
    return ControlPoints

def make_downstream_points(fn):
    DownstreamPoints = []
    with open(fn,'r') as f:
        next(f)
        for line in f:
            line = line.strip('\n').split(',')
            DownstreamPoints.append(DownstreamPoint(line))

    return DownstreamPoints




def orifice_xsect_grab(controlDict,orifices):
    # Add items from orifices dictionary to the control dict
    for i in controlDict:
        try:
            controlDict[i].update(orifices[i])
        except:
            pass
        
def pump_curve_grab(controlDict, pumps):
    # Add information to controlDict with pump curve info.
    for i in controlDict:
        try:
            controlDict[i].update(pumps[i])
        except:
            pass

def get_depth(elements,conduitDict,storageDict):
    for element in elements:
        if elements[element]['type'] == 'link':
            elements[element]['max_depth'] = conduitDict[element]['geom1']
        elif elements[element]['type'] == 'storage':
            elements[element]['max_depth'] = storageDict[element]['max_depth']
        else:
            pass

def get_q_full_and_other(elements,conduits,storages,junctions):
    for element in elements:
        if elements[element]['type'] == 'link':
            elements[element]['max_flow'] = conduits[element]['q_full']
        elif elements[element]['type'] == 'storage':
            elements[element]['total_storage'] = storages[element]['total_storage']
        elif elements[element]['type'] == 'junction':
            elements[element]['max_depth'] = junctions[element]['max_depth']
        else:
            pass

# SAVING

def push_meta(run,outF,ControlPoints,DownstreamPoints,derivative):
    # Make string to push to csv metadata file
    push_list = []

    push_list.append(run.flood_count)
    push_list.append(outF)
    for g in range(1,run.groups+1):
        print(g+1)
        push_list.append('Group' + str(g))
        
        push_list.append('UP')
        cps = [c for c in ControlPoints if c.group == g]
        for c in cps:
            push_list.append(c.u_name)
            push_list.append(c.u_param)
            
            if derivative:
                push_list.append(c.ds_param)
        
        push_list.append('DOWN')
        dps = [d for d in DownstreamPoints if d.group == g]
        for d in dps:
            push_list.append(d.d_name)
            push_list.append(d.epsilon)
            push_list.append(d.set_point)
            
            if derivative:
                push_list.append(d.gamma)
                push_list.append(d.set_derivative)

    push_str = ','.join(map(str, push_list))
    return push_str


def time_nano(dt):
    epoch0 = datetime.datetime(1970,1,1)
    epoch0 = epoch0.replace(tzinfo = pytz.UTC)

    return str(int((dt.replace(tzinfo = pytz.UTC) - epoch0).total_seconds()*1000000000))   

# PLOTTTING Stuff
def make_extract_string(name,el_type,measure):

    node_keys = {'depth': 'Depth_above_invert',
                'head': 'Hydraulic_head',
                'volume': 'Volume_stored_ponded',
                'inflow_lat': 'Lateral_inflow',
                'flow': 'Total_inflow',
                'flooding': 'Flow_lost_flooding'}
    link_keys = {'flow': 'Flow_rate',
                'depth': 'Flow_depth',
                'velocity': 'Flow_velocity',
                'froude': 'Froude_number',
                'cap': 'Capacity'}
        

    if el_type == 'storage' or el_type == 'junction':
        el_type = 'node'
        measure = node_keys[measure]
    else:
        measure = link_keys[measure]
        #el_type already link
        
    return el_type + ',' + name + ',' + measure

