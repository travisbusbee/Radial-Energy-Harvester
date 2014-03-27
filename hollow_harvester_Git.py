from mecode import G
import math
import os


import numpy as np

# Travis' Computer Outfile
#outfile = r"C:\Users\tbusbee\Documents\GitHub\Radial-Energy-Harvester\radial-harvester-testing.txt"

# Robomama Outfile
outfile = r"C:\Users\Lewis Group\Documents\GitHub\Radial-Energy-Harvester\radial-harvester-testing.pgm"

#calfile =  r"C:\Users\Lewis Group\Desktop\Busbee\profilometer_output_030214_1.txt"

alignment_file_1 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_1.txt"
alignment_file_2 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_2.txt"
alignment_file_3 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_3.txt"
alignment_file_4 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_4.txt"

cal_data = None#load_and_curate(calfile, reset_start=(2, -2))

g = G(
    outfile=outfile,
    #header=None,
    #footer=None,
    #cal_data=cal_data,
    print_lines=False,
    )
    
zero=(79.216500, 82.533100, 60, 93.852950)
pressure_box = 4
g.setup()
#z_start = 1.6
def pressure_purge():
    g.toggle_pressure(pressure_box)
    g.write('$DO6.0=1')
    g.dwell(0.75)
    g.write('$DO6.0=0')
    g.toggle_pressure(pressure_box)
    g.dwell(0.5)

def spiral(radius, over, x_center, y_center,  direction = 'CW'):
    g.abs_move(x=x_center, y=(y_center-over*0.5))
    g.set_valve(num = 0, value = 1)
    repeats = 2*(radius/over)
    count = 0
    
    negative=-1
    for i in range(int(repeats)):
        count=count+1
        power = count + 1
        sign=pow(negative, power)
        
    
        if sign >0:
            g.arc(direction=direction, radius = 0.5*(over*count), x=0, y=(over*count))
        else:
            g.arc(direction=direction, radius = 0.5*(over*count), x=0, y=(-over*count))
    g.set_valve(num = 0, value = 0)

def hollow_spiral(inner_radius, outer_radius, over, x_center, y_center, direction = 'CW'):
    r_in = inner_radius + 0.5*over
    r_out = outer_radius - 0.5*over
    r_dif = r_out-r_in
    repeats = 2*(r_dif/over)
    g.set_valve(num = 0, value = 0)
    g.abs_move(x=x_center, y=(y_center-r_in))
    g.set_valve(num = 0, value = 1)
    count = 0
    negative=-1
    for i in range(int(repeats)):
        count=count+1
        power = count + 1
        sign=pow(negative, power)
    
        if sign >0:
            g.arc(direction=direction, radius = (inner_radius + 0.5*(over*count)), x=0, y=(2*inner_radius+(over*count)))
        else:
            g.arc(direction=direction, radius = (inner_radius + 0.5*(over*count)), x=0, y=-(2*inner_radius+(over*count)))
    g.set_valve(num = 0, value = 0)
    
def stacked_hollow_spirals(min_radius, max_radius, outer_radius, r_step, layer_height, over, x_center, y_center, direction='CW', nozzle = 'A', taper = 'min_to_max'):
    layers=(max_radius-min_radius)/r_step
    upper = (2*np.pi - (((0.5*np.pi)/layers)/2))
    lower = 1.5*np.pi + (((0.5*np.pi)/layers)/2)
    angle_step = (0.5*np.pi)/layers
    if taper=='min_to_max':
        angle=np.arange(lower, upper, angle_step)
    elif taper == 'max_to_min':
        angle=np.arange(upper, lower, angle_step) 
    else:
        raise RuntimeError('Why dont you make like a tree and go fuck yourself: {}'.format(taper))
    
    R_in = np.cos(angle)*(max_radius-min_radius)+min_radius
    
   
    g.write('; There are {} stacked hollow layers'.format(layers))
    for i in range(int(layers)):
        hollow_spiral(inner_radius=R_in[i], outer_radius = outer_radius, over=over, x_center=x_center, y_center=y_center, direction = direction)
        g.move(**{nozzle:layer_height})
    

def recalculate_xy_intercepts(x_center, y_center, over, center_to_top, xr, xl, temporary_radius, counter):
    temp_inner = temporary_radius + counter*over               
    x_right = xr + counter*over
    x_left = xl - counter*over
    x_over = x_right - x_center
    y_start = y_center + center_to_top 
    theta = math.acos(x_over/temp_inner)
    y_end = math.tan(theta)*x_over + y_center         


def arc_meander(x_center, y_center, R_min, R_max, over, start_angle, arc_length, start_direction = 'CW', start = 'R_min'):
    start_angle = math.radians(start_angle)
    arc_length = math.radians(arc_length)
    temp_inner = R_min + 0.5*over
    temp_outer = R_max - 0.5*over
    loops = ((temp_outer - temp_inner)/over)/2
    print loops
    if start=='R_min':
        R_start = temp_inner
        sign = 1
    elif start=='R_max':
        R_start = temp_outer
        sign = -1
    else:
        raise RuntimeError('invalid starting radius: {}'.format(start))
    x_start = x_center + math.cos(start_angle)*R_start
    y_start = y_center + math.sin(start_angle)*R_start
    count = 0
    
    if start_direction == 'CW':
        turnaround_angle = start_angle - arc_length
        direction1 = 'CW'
        direction2 = 'CCW'
    elif start_direction == 'CCW':
        turnaround_angle = start_angle + arc_length
        direction1 = 'CCW'
        direction2 = 'CW'
    else:
        raise RuntimeError('invalid start direction: {}'.format(start_direction))
    
    
    for i in range(int(loops)):
        R_working = R_start + count*over*sign
        x1 = x_center + math.cos(start_angle)*R_working
        y1 = y_center + math.sin(start_angle)*R_working
        x2 = x_center + math.cos(turnaround_angle)*R_working
        y2 = y_center + math.sin(turnaround_angle)*R_working
        if count > 0:
            g.abs_move(x=x1, y=y1)
        
        g.abs_move(x1, y1)
        g.abs_arc(direction=direction1, radius = R_working, x=x2, y=y2)
        count = count + 1
        R_working = R_start + count*over*sign
        x2 = x_center + math.cos(turnaround_angle)*R_working
        y2 = y_center + math.sin(turnaround_angle)*R_working
        x1 = x_center + math.cos(start_angle)*R_working
        y1 = y_center + math.sin(start_angle)*R_working
        g.abs_move(x2, y2)
        g.abs_arc(direction = direction2, radius = R_working, x=x1, y=y1)
        count = count + 1
    
def print_insulating_meanders():    
   arc_meander(x_center, y_center, R_min, R_max, over, start_angle, arc_length, start_direction = 'CW', start = 'R_min')
   
   
   
def print_tailed_ring(over, xr_in, xl_in, y_start, y_end_in, r_in, r_out, xr_out, xl_out, y_end_out, layers, layer_height, speed, pressure, nozzle = 'A'): 
    
    for i in range(int(layers)):
        g.abs_move(x=xr_in, y=y_start)
        g.abs_move(x=xr_in, y=y_end_in)
        g.abs_arc(direction='CW', radius = -r_in, x=xl_in, y=y_end_in)
        g.abs_move(x=xl_in, y=y_start)
        g.abs_move(x=xl_out)                     
        g.abs_move(x=xl_out, y=y_end_out)
        g.abs_arc(direction='CCW', radius = -r_out, x=xr_out, y=y_end_out)
        g.abs_move(x=xr_out, y=y_start)
        g.abs_move(x=xr_in)
        g.move(**{nozzle:layer_height})
    
                                                                                                                                    
def inner_electrode(x_center, y_center, x_offset, inner_radius, outer_radius, over, center_to_top, speed, layers, layer_height, pressure, nozzle ):
    count = 1
    temp_inner = inner_radius + count*0.5*over               
    x_right = x_center + count*0.5*over+ x_offset
    x_left = x_center - count*0.5*over - x_offset
    x_over = x_right - x_center
    y_start = y_center + center_to_top 
    theta = math.acos(x_over/temp_inner)
    y_end = math.tan(theta)*x_over + y_center
    temp_inner_outer = temp_inner + count*over               
    x_right_outer = x_right + count*over
    x_left_outer = x_left - count*over
    x_over_outer = x_right_outer - x_center
    theta = math.acos(x_over_outer/temp_inner_outer)
    y_end_outer = math.tan(theta)*x_over_outer + y_center
    g.feed(speed)
    print_tailed_ring(over=over, xr_in=x_right, xl_in=x_left, y_start=y_start, y_end_in=y_end, r_in=temp_inner, r_out=temp_inner_outer, xr_out=x_right_outer, 
                        xl_out=x_left_outer, y_end_out=y_end_outer, layers = layers, layer_height = layer_height, speed = speed, pressure = pressure, nozzle = nozzle)
   
def single_straight_stack(x, layer_height, z, layers, speed, pressure, nozzle, valve):
   
    count = 0
    for i in range(layers):
        count = count +1
        g.set_valve(num = valve, value = 1)
        g.move(x=x)
        g.move(**{nozzle:layer_height})
        g.move(x=-x)
        if count != layers:
            g.move(**{nozzle:layer_height})
        else:
            g.set_valve(num = valve, value = 0)
            g.clip(axis=nozzle, direction='-x', height=4)
    
def double_straight_stack(x, over, layer_height, z, layers, speed, pressure, nozzle, valve):
   
    count = 0
    for i in range(layers):
        count = count +1
        g.set_valve(num = valve, value = 1)
        g.move(x=x)
        g.move(y=over)
        g.move(x=-x)
        g.move(y=-over)
        if count != layers:
            g.move(**{nozzle:layer_height})
        else:
            g.set_valve(num = valve, value = 0)
            g.dwell(0.2)
            g.clip(axis=nozzle, direction='-x', height=2)
        
def print_straight_electrode(x, z, tail, layer_height, layers, speed, pressure, nozzle, valve):
    g.feed(15)
    g.set_pressure(com_port=pressure_box, value=pressure)    
    g.abs_move(**{nozzle:z})
    g.feed(speed)
    g.set_valve(num = valve, value = 1)
    g.dwell(0.2)
    g.move(x=tail)
    single_straight_stack(x, layer_height, z, layers, speed, pressure, nozzle, valve)  
    g.set_valve(num = valve, value = 0)

def print_double_electrode(x, over, z, tail, layer_height, layers, speed, pressure, nozzle, valve):
    g.feed(15)
    g.set_pressure(com_port=pressure_box, value=pressure)    
    g.abs_move(**{nozzle:z})
    g.feed(speed)
    g.set_valve(num = valve, value = 1)
    g.dwell(0.8)
    #g.move(x=tail)
    double_straight_stack(x, over, layer_height, z, layers, speed, pressure, nozzle, valve)  
    g.set_valve(num = valve, value = 0)      
                  
def partial_tailed_ring(x_center, y_center, z_start, x_offset, inner_radius, over, center_to_top, layers, layer_height, speed, pressure, orientation = 'tail_up', side = 'left', angle = 25, nozzle = 'A'):
    temp_inner = inner_radius + 0.5*over  
    temp_outer = inner_radius + 1.5*over  
    x_over_in = x_offset + 0.5*over 
    x_over_out = x_over_in + over
    radians = math.radians(angle)
    x_turn_in = temp_inner*math.sin(radians)
    x_turn_out = temp_outer*math.sin(radians)
    y_turn_in = temp_inner*math.cos(radians)
    y_turn_out = temp_outer*math.cos(radians)              
    if side == 'left':
        x_in = x_center - x_offset  - 0.5*over
        x_out = x_in - over
        x_turn_in = x_center - x_turn_in
        x_turn_out = x_center - x_turn_out
        indicator = 1
    elif side == 'right':
        x_in = x_center + x_offset + 0.5*over
        x_out = x_in + over
        x_turn_in = x_center + x_turn_in
        x_turn_out = x_center + x_turn_out
        indicator = -1
    else:
        raise RuntimeError('invalid nside: {}'.format(side))          

    theta_in = math.acos(x_over_in/temp_inner)
    theta_out = math.acos(x_over_out/temp_outer)
    y_dist_in = math.tan(theta_in)*x_over_in
    y_dist_out = math.tan(theta_out)*x_over_out 
    direction1 = 'CW'
    direction2 = 'CCW'                                 
    if orientation == 'tail_up': 
        y_start = y_center + center_to_top
        y_end_in = y_center + y_dist_in
        y_end_out = y_center + y_dist_out
        y_turn_in = y_center - y_turn_in
        y_turn_out = y_center - y_turn_out
        if indicator > 0:
            direction1 = 'CCW'
            direction2 = 'CW'
    else:
        y_start = y_center - center_to_top                     
        y_end_in = y_center - y_dist_in
        y_end_out = y_center - y_dist_out                    
        y_turn_in = y_center + y_turn_in
        y_turn_out = y_center + y_turn_out                        
        if indicator < 0:
           direction1 = 'CCW'
           direction2 = 'CW' 
    g.feed(speed)  
    g.abs_move(x=x_in, y=y_start)
    g.abs_move(**{nozzle:z_start})
    for i in range(int(layers)):
        g.abs_move(x=x_in, y=y_start)
        g.set_valve(num = 0, value = 1)
        g.abs_move(x=x_in, y=y_end_in)
        g.abs_arc(direction=direction1 , radius = temp_inner, x=x_turn_in, y=y_turn_in) 
        g.abs_move(x=x_turn_out, y=y_turn_out) 
        g.abs_arc(direction = direction2 , radius = temp_outer, x=x_out, y= y_end_out)                               
        g.abs_move(x=x_out, y=y_start)
        g.abs_move(x=x_in, y=y_start)
        g.move(**{nozzle:layer_height})                                   
    g.set_valve(num = 0, value = 0)
def nozzle_change_vars(nozzles = 'ab'):
    g.feed(40)
    g.home()
    if g.cal_data == None:
        cal_off = True
    else:
        cal_off = False
    g.cal_data=None
    g.dwell(0.25)
    g.write(';----------nozzle change------------')
    if nozzles=='ab':
        g.abs_move(A=50)
        g.write('G1 X($Bx-$Ax-($Bx_dif-$Ax_dif))  Y($Ay-$By+($Ay_dif-$By_dif))')
    elif nozzles=='ac':
        g.abs_move(A=50)
        g.write('G1 X($Cx-$Ax-($Cx_dif-$Ax_dif))  Y($Ay-$Cy+($Ay_dif-$Cy_dif))')    
    elif nozzles=='ad':
        g.abs_move(A=50)
        g.write('G1 X($Dx-$Ax-($Dx_dif-$Ax_dif))  Y($Ay-$Dy+($Ay_dif-$Dy_dif))')
    elif nozzles=='ba':
        g.abs_move(B=50)
        g.write('G1 X($Ax-$Bx-($Ax_dif-$Bx_dif))  Y($By-$Ay+($By_dif-$Ay_dif))')  
    elif nozzles=='bc':
        g.abs_move(B=50)
        g.write('G1 X($Cx-$Bx-($Cx_dif-$Bx_dif))  Y($By-$Cy+($By_dif-$Cy_dif))')
    elif nozzles=='bd':
        g.abs_move(B=50)
        g.write('G1 X($Dx-$Bx-($Dx_dif-$Bx_dif))  Y($By-$Dy+($By_dif-$Dy_dif))')
    elif nozzles=='ca':
        g.abs_move(C=50)
        g.write('G1 X($Ax-$Cx-($Ax_dif-$Cx_dif))  Y($Cy-$Ay+($Cy_dif-$Ay_dif))')
    elif nozzles=='cb':
        g.abs_move(C=50)
        g.write('G1 X($Bx-$Cx-($Bx_dif-$Cx_dif))  Y($Cy-$By+($Cy_dif-$By_dif))')
    elif nozzles=='cd':
        g.abs_move(C=50)
        g.write('G1 X($Dx-$Cx-($Dx_dif-$Cx_dif))  Y($Cy-$Dy+($Cy_dif-$Dy_dif))')
    elif nozzles=='da':
        g.abs_move(D=50)
        g.write('G1 X($Ax-$Dx-($Ax_dif-$Dx_dif))  Y($Dy-$Ay+($Dy_dif-$Ay_dif))')
    elif nozzles=='db':
        g.abs_move(D=50)
        g.write('G1 X($Bx-$Dx-($Bx_dif-$Dx_dif))  Y($Dy-$By+($Dy_dif-$By_dif))')
    elif nozzles=='dc':
        g.abs_move(D=50)
        g.write('G1 X($Cx-$Dx-($Cx_dif-$Dx_dif))  Y($Dy-$Cy+($Dy_dif-$Cy_dif))')
    else:
        g.write('; ---------- input a real nozzle change input...ya idiot--------')   
    if cal_off == False:
        g.cal_data=load_and_curate(calfile, reset_start=(2, -2))       
        g.cal_axis = nozzles[1].upper()         

def calculate_relative_z(reference_nozzle = 'A'):
    if reference_nozzle == 'A':
        g.write('$zA = -{}' .format(zero[0]))
        g.write('$zB = $zA + ($zMeasureB - $zMeasureA) + ($Az_dif-$Bz_dif)')
        g.write('$zC = $zA + ($zMeasureC - $zMeasureA) + ($Az_dif-$Cz_dif)')
        g.write('$zD = $zA + ($zMeasureD - $zMeasureA) + ($Az_dif-$Dz_dif)')
    elif reference_nozzle == 'B':
            g.write('$zB = -{}' .format(zero[1]))
            g.write('$zA = $zB + ($zMeasureA - $zMeasureB) + ($Bz_dif-$Az_dif)')
            g.write('$zC = $zB + ($zMeasureC - $zMeasureB) + ($Bz_dif-$Cz_dif)')
            g.write('$zD = $zB + ($zMeasureD - $zMeasureB) + ($Bz_dif-$Dz_dif)')
    elif reference_nozzle == 'C':
            g.write('$zC = -{}' .format(zero[2]))
            g.write('$zA = $zC + ($zMeasureA - $zMeasureC) + ($Cz_dif-$Az_dif)')
            g.write('$zB = $zC + ($zMeasureB - $zMeasureC) + ($Cz_dif-$Bz_dif)')
            g.write('$zD = $zC + ($zMeasureD - $zMeasureC) + ($Cz_dif-$Dz_dif)')
    elif reference_nozzle == 'D':
            g.write('$zD = -{}' .format(zero[3]))
            g.write('$zA = $zD + ($zMeasureA - $zMeasureD) + ($Dz_dif-$Az_dif)')
            g.write('$zB = $zD + ($zMeasureB - $zMeasureD) + ($Dz_dif-$Bz_dif)')
            g.write('$zC = $zD + ($zMeasureC - $zMeasureD) + ($Dz_dif-$Cz_dif)')

def set_home_in_aerotech():
    g.write('G92 A(-$zA-5) B(-$zB-5) C(-$zC - 5) D(-$zD - 5)')   
    
def recall_alignment(nozzle = 'A'):
   if nozzle=='A':
        g.write(open(alignment_file_1).read()) 
   elif nozzle=='B':
        g.write(open(alignment_file_2).read()) 
   elif nozzle=='C':
        g.write(open(alignment_file_3).read())
   elif nozzle=='D':
        g.write(open(alignment_file_4).read())
   elif nozzle =='all':
        g.write(open(alignment_file_1).read())
        g.write(open(alignment_file_2).read())
        g.write(open(alignment_file_3).read())
        g.write(open(alignment_file_4).read())

def print_top_inner_electrode(z_start, nozzle):
    g.abs_move(**{nozzle:z_start})
    g.set_valve(num = 0, value = 1)
    inner_electrode(x_center = 0, y_center=0, x_offset=0, inner_radius = 7.5, outer_radius = 20, over=0.5, center_to_top = 20, speed = 8, layers = 5, layer_height = 0.3, pressure = 40, nozzle = nozzle)  
    g.set_valve(num = 0, value = 0)                                                                                                  
def print_top_outer_electrode(z_start, nozzle):
    
    partial_tailed_ring(x_center=0, y_center=0, z_start=z_start, x_offset= 0.75, inner_radius = 12.5, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_up', side = 'left', angle = 25, nozzle = nozzle)  
    g.move(**{nozzle:5})                           
    partial_tailed_ring(x_center=0, y_center=0, z_start=z_start, x_offset= 0.75, inner_radius = 12.5, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_up', side = 'right', angle = 25, nozzle = nozzle)
    g.move(**{nozzle:5})
def print_bottom_inner_electrode(z_start, nozzle):
    
    partial_tailed_ring(x_center=0, y_center=0, z_start=z_start, x_offset= 0, inner_radius = 10, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'left', angle = 25, nozzle = nozzle)
    g.move(**{nozzle:5})
    partial_tailed_ring(x_center=0, y_center=0, z_start=z_start, x_offset= 0, inner_radius = 10, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'right', angle = 25, nozzle = nozzle)
    g.move(**{nozzle:5})

def print_bottom_outer_electrode(z_start, nozzle):
    
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.75, z_start=z_start, inner_radius = 15, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'left', angle = 25, nozzle = nozzle)
    g.move(**{nozzle:5})
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.75, z_start=z_start, inner_radius = 15, over = 0.5, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'right', angle = 25, nozzle = nozzle)
    g.move(**{nozzle:5})



recall_alignment(nozzle = 'all')
g.feed(30)
#g.align_zero_nozzle(nozzle='A', floor=-49.25, deltafast=0.85, deltaslow=0.1, start=-15)
#g.align_zero_nozzle(nozzle='B', floor=-49.25, deltafast=0.85, deltaslow=0.1, start=-15)
#g.align_zero_nozzle(nozzle='D', floor=-49.25, deltafast=0.85, deltaslow=0.1, start=-15)
g.save_alignment(nozzle = 'all')
g.feed(30)
g.abs_move(A=-2, B=-2, C=-2, D=-2)
g.abs_move(x=352, y=70.24)#197.96
g.write('G1 X$Ax_dif  Y$Ay_dif')
g.set_home(x=0, y=0)
pressure_purge()
g.toggle_pressure(pressure_box)
calculate_relative_z(reference_nozzle = 'B')

g.abs_move(A=-5, B=-5, C=-5, D=-5)
set_home_in_aerotech()

nozzle_change_vars('ab')
g.set_home(x=0, y=0)

#print_straight_electrode(x = 20, z = 0.15, tail = 10, layer_height = 0.15, layers = 5, speed = 4, pressure = 45, nozzle = 'B', valve = 1)
g.feed(15)

g.abs_move(x=20, y=-20)
g.set_pressure(pressure_box, 50)
g.abs_move(B=0.2)
g.feed(6)
g.set_valve(num = 1, value = 1)
g.meander(x=20, y=-20, spacing=0.34, orientation = 'y', tail = False)
g.set_valve(num = 1, value = 0)
#print_double_electrode(x = 20, over=-0.34, z = 0.22, tail = 10, layer_height = 0.25, layers = 8, speed = 6, pressure = 40, nozzle = 'B', valve = 1)

#g.move(y=-2.2)
#print_double_electrode(x = 20, over=-0.34, z = 0.22, tail = 10, layer_height = 0.25, layers = 8, speed = 6, pressure = 40, nozzle = 'B', valve = 1)

#hollow_spiral(inner_radius = 5, outer_radius = 20, over = 1 , x_center = 0, y_center = 0, direction = 'CW')

g.toggle_pressure(pressure_box)
#
#arc_meander(x_center=0, y_center=0, R_min = 8.9, R_max = 10.85, over = 0.65, start_angle = 80, arc_length = 12, start_direction = 'CW', start = 'R_min')
#arc_meander(x_center=0, y_center=0, R_min = 10.85, R_max = 12.5, over = 0.65, start_angle = 83, arc_length = 15, start_direction = 'CW', start = 'R_min')
#g.abs_move(z=0.2)

#g.abs_move(z=0.4)
#g.set_pressure(pressure_box, 19)
#g.feed(12)
#pressure_purge()
#g.abs_move(A=0.2)
#spiral(radius=10, over=0.5, x_center=0, y_center=0, direction = 'CW')
#g.abs_move(A=0.5)
#g.set_valve(num = 0, value = 1)
#stacked_hollow_spirals(min_radius = 2, max_radius=4, outer_radius = 10, r_step=0.25, layer_height = 0.35, over=0.5, x_center=0, y_center=0, direction='CW', nozzle = 'A')
#g.set_valve(num = 0, value = 0)




#print_top_inner_electrode(z_start=0.2, nozzle = 'A')
#print_top_outer_electrode(z_start=0.2, nozzle = 'A')
#print_bottom_inner_electrode(z_start=0.2, nozzle = 'A')
#print_bottom_outer_electrode(z_start=0.2, nozzle = 'A')
#
#
#
#g.move(A=40)












        
#spiral(radius = 30, over=0.5, x_center = 0, y_center = 0, direction = 'CW')


g.teardown()