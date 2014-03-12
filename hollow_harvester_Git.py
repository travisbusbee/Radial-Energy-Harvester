from mecode import G
import math
import os


import numpy as np

# Travis' Computer Outfile
outfile = r"C:\Users\tbusbee\Documents\GitHub\Radial-Energy-Harvester\radial-harvester-testing.txt"
#calfile =  r"C:\Users\Lewis Group\Desktop\Busbee\profilometer_output_030214_1.txt"

#alignment_file_1 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_1.txt"
#alignment_file_2 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_2.txt"
#alignment_file_3 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_3.txt"
#alignment_file_4 = r"C:\Users\Lewis Group\Desktop\Alignment\alignment_values_4.txt"

cal_data = None#load_and_curate(calfile, reset_start=(2, -2))

g = G(
    outfile=outfile,
    #header=None,
    #footer=None,
    #cal_data=cal_data,
    print_lines=False,
    )
g.setup()

def spiral(radius, over, x_center, y_center, direction = 'CW'):
    g.abs_move(x=x_center, y=(y_center-over*0.5))
    repeats = 2*(radius/over)
    count = 0
    
    negative=-1
    for i in range(int(repeats)):
        count=count+1
        power = count + 1
        sign=pow(negative, power)
        
        g.write('sign is {}'.format(sign))
        if sign >0:
            g.arc(direction=direction, radius = 0.5*(over*count), x=0, y=(over*count))
        else:
            g.arc(direction=direction, radius = 0.5*(over*count), x=0, y=(-over*count))


def hollow_spiral(inner_radius, outer_radius, over, x_center, y_center, direction = 'CW'):
    r_in = inner_radius + 0.5*over
    r_out = outer_radius - 0.5*over
    r_dif = r_out-r_in
    repeats = 2*(r_dif/over)
    g.abs_move(x=x_center, y=(y_center-r_in))
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

def stacked_hollow_spirals(min_radius, max_radius, outer_radius, r_step, layer_height, over, x_center, y_center, direction='CW', nozzle = 'A'):
    layers=(max_radius-min_radius)/r_step
    upper = (2*np.pi - (((0.5*np.pi)/layers)/2))
    lower = 1.5*np.pi + (((0.5*np.pi)/layers)/2)
    angle_step = (0.5*np.pi)/layers
    angle=np.arange(lower, upper, angle_step)
    
    R_in = np.cos(angle)*(max_radius-min_radius)+min_radius
    
   
    
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
    
                                                                                                                                    
def inner_electrode(x_center, y_center, x_offset, inner_radius, outer_radius, over, center_to_top, speed, layers, layer_height, pressure ):
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
                        xl_out=x_left_outer, y_end_out=y_end_outer, layers = layers, layer_height = layer_height, speed = speed, pressure = pressure, nozzle = 'z')
   

def partial_tailed_ring(x_center, y_center, x_offset, inner_radius, over, center_to_top, layers, layer_height, speed, pressure, orientation = 'tail_up', side = 'left', angle = 25, nozzle = 'A'):
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
    for i in range(int(layers)):
        g.abs_move(x=x_in, y=y_start)
        g.abs_move(x=x_in, y=y_end_in)
        g.abs_arc(direction=direction1 , radius = temp_inner, x=x_turn_in, y=y_turn_in) 
        g.abs_move(x=x_turn_out, y=y_turn_out) 
        g.abs_arc(direction = direction2 , radius = temp_outer, x=x_out, y= y_end_out)                               
        g.abs_move(x=x_out, y=y_start)
        g.abs_move(x=x_in, y=y_start)
        g.move(**{nozzle:layer_height})                                   

def print_top_inner_electrode():
    g.abs_move(z=0.15)
    inner_electrode(x_center = 0, y_center=0, x_offset=0, inner_radius = 10, outer_radius = 20, over=0.6, center_to_top = 20, speed = 8, layers = 5, layer_height = 0.3, pressure = 40)  
                                                                                                      
def print_top_outer_electrode():
    g.abs_move(z=0.15)
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.9, inner_radius = 13.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_up', side = 'left', angle = 25, nozzle = 'z')  
    g.abs_move(z=0.15)                            
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.9, inner_radius = 13.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_up', side = 'right', angle = 25, nozzle = 'z')

def print_bottom_inner_electrode():
    g.abs_move(z=0.15)
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0, inner_radius = 12.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'left', angle = 25, nozzle = 'z')
    g.abs_move(z=0.15)
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0, inner_radius = 12.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'right', angle = 25, nozzle = 'z')


def print_bottom_outer_electrode():
    g.abs_move(z=0.15)
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.9, inner_radius = 15.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'left', angle = 25, nozzle = 'z')
    g.abs_move(z=0.15)
    partial_tailed_ring(x_center=0, y_center=0, x_offset= 0.9, inner_radius = 15.5, over = 0.6, center_to_top = 20, layers = 5, layer_height = 0.3, speed = 8, pressure = 40, orientation = 'tail_down', side = 'right', angle = 25, nozzle = 'z')



#hollow_spiral(inner_radius = 5, outer_radius = 20, over = 1 , x_center = 0, y_center = 0, direction = 'CW')



arc_meander(x_center=0, y_center=0, R_min = 11, R_max = 40, over = 0.8, start_angle = 88.5, arc_length = 22, start_direction = 'CW', start = 'R_min')


#stacked_hollow_spirals(min_radius = 2, max_radius=10, outer_radius = 25, r_step=1, layer_height = 1, over=1, x_center=0, y_center=0, direction='CW', nozzle = 'z')
#print_top_inner_electrode()
#print_top_outer_electrode()
#print_bottom_inner_electrode()
#print_bottom_outer_electrode()
















        
#spiral(radius = 30, over=0.5, x_center = 0, y_center = 0, direction = 'CW')


g.teardown()