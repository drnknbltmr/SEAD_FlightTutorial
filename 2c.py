import numpy as np

#airplane info
b = 26.17 
aircraft_length = 39.13
fuselage_length = 36.466 #m

EOW = 23178 #kg
Wing_w = 0.197*EOW
HTail_w = 0.025*EOW
VTail_w = 0.018*EOW
Fuselage_w = 0.35*EOW
MLG_w = 0.058 * EOW
NLG_w = 0.008*EOW
Propulsion_sys_w = 0.133 * EOW
Cockpit_sys_w = 0.023*EOW

Wing_group_W = Wing_w + MLG_w
Fuselage_group_W = Fuselage_w + NLG_w + Cockpit_sys_w + Propulsion_sys_w + HTail_w + VTail_w
Other_comp_w =EOW - (Wing_group_W + Fuselage_group_W)



#components positions - change them after from front
#Wing cg position
c_t = 1.31
c_r = 5.69
x_to_c_t = (4.2 * 247.16)/100 
fuselage_diam = 2.695
c_035b = 0.65 * c_r + 0.35 * c_t
front_spar = 0.2* c_035b # estimate from adsee 1
rear_spar = 0.7 * c_035b # estimate from adsee 1
x_cg_c = front_spar + 0.7 * (rear_spar - front_spar)
front_sweep = 30 #degree
Wing_position_front = 39.133 - 5.11 - 16.97 - fuselage_diam * np.tan(front_sweep/180*np.pi) #m from nose
x_cg_wing = Wing_position_front + x_cg_c + (0.35*b/2) * np.tan(front_sweep/180*np.pi)


#Horizontal tail cg position
b_ht = 8.541
front_sweep_ht = 34 #degree
back_sweep_ht = 17 #degree
length = 4.04
c_r_ht = 4.04 - b_ht/2 * np.tan(back_sweep_ht/180*np.pi)
c_t_ht = (0.5 * 247.16)/100
c038_ht = 0.55 * c_t_ht + 0.45 * c_r_ht
hor_tail_position_front = 39.133 - 4.04

c_cg_ht = 0.42 * c038_ht
x_cg_HTail = hor_tail_position_front + c_cg_ht + (0.55*b_ht/2) * np.tan(front_sweep_ht/180*np.pi)

#vertical tail cg position
l_vt = 1.4*2
front_sweep_vt = 39 #degree
back_sweep_vt = 30 #degree
c_r_vt = 1.3 * 2
c_t_vt = 0.9 * 2
c038_vt = 0.55 * c_t_vt + 0.45 * c_r_vt
ver_tail_position_front = 28.6+1.09+1.2*2
c_cg_vt = 0.42 * c038_vt
x_cg_Vertical_tail = ver_tail_position_front + c_cg_vt + (0.55*l_vt) * np.tan(front_sweep_vt/180*np.pi)

#fuselage cg position
x_cg_fuselage = 0.47*fuselage_length
print("fuselage cg from nose: ", x_cg_fuselage)
# landing gears cg position

x_cg_MLG = aircraft_length-16.97-0.1*5.11 #based on rare wing spar distance

x_cg_NLG = 1.1*2

#propiulsion system and cockpit system cg position (nacelle = engine cg)
x_cg_propulsion = 0.4*1.33*2.4716 + 28.6 

#cockpit system cg position
x_cg_cockpit = x_cg_NLG 

Wing_grpoup_cg = (Wing_w * x_cg_wing + MLG_w * x_cg_MLG) / (Wing_w + MLG_w)
print("Wing group cg from nose: ", Wing_grpoup_cg)
Fuselage_group_cg = (Fuselage_w * x_cg_fuselage + NLG_w * x_cg_NLG + Cockpit_sys_w * x_cg_cockpit + Propulsion_sys_w * x_cg_propulsion + HTail_w * x_cg_HTail + VTail_w * x_cg_Vertical_tail) / Fuselage_group_W
print("Fuselage group cg from nose: ", Fuselage_group_cg)
EOW_cg = (Wing_group_W * Wing_grpoup_cg + Fuselage_group_W * Fuselage_group_cg) / (Wing_group_W + Fuselage_group_W)

c_root = 5.7
taper = c_t/c_root
#x_MAC = (2*c_root)/3 * (1+taper+taper**2)/(1+taper)
x_MAC = 3.48 #m
y_mac = (x_to_c_t + fuselage_diam/2)/3 * (1+2*taper)/(1+taper)
#x_LEMAC = Wing_position_front + y_mac * np.tan(front_sweep/180*np.pi)
x_LEMAC = 22.866 - 3.6576
#components position based on LEMAC
x_cg_LEMAC_wing = x_cg_wing - x_LEMAC
x_cg_LEMAC_HTail = x_cg_HTail - x_LEMAC
x_cg_LEMAC_Vertical_tail = x_cg_Vertical_tail - x_LEMAC
x_cg_LEMAC_fuselage = x_cg_fuselage - x_LEMAC
x_cg_LEMAC_MLG = x_cg_MLG - x_LEMAC
x_cg_LEMAC_NLG = x_cg_NLG - x_LEMAC
x_cg_LEMAC_propulsion = x_cg_propulsion - x_LEMAC
x_cg_LEMAC_cockpit = x_cg_cockpit - x_LEMAC

Wing_grpoup_cg_LEMAC = (Wing_w * x_cg_LEMAC_wing + MLG_w * x_cg_LEMAC_MLG) / Wing_group_W
Fuselage_group_cg_LEMAC = (Fuselage_w * x_cg_LEMAC_fuselage + NLG_w * x_cg_LEMAC_NLG + Cockpit_sys_w * x_cg_LEMAC_cockpit + Propulsion_sys_w * x_cg_LEMAC_propulsion + HTail_w * x_cg_LEMAC_HTail + VTail_w * x_cg_LEMAC_Vertical_tail) / Fuselage_group_W
EOW_cg_LEMAC = (Wing_group_W * Wing_grpoup_cg_LEMAC + Fuselage_group_W * Fuselage_group_cg_LEMAC) / (Wing_group_W + Fuselage_group_W)   
print("EOW cg from nose: ", EOW_cg)