import pandas as pd
import numpy as np
import math

# Function to read from dumpfile and create dataFrame
def ReadDumpFile(filename):
    filepath = '../' + filename
    col_names = ['id', 'type', 'diameter', 'mass', 'xu', 'yu', 'zu','fx', 'fy', 'fz']
    df = pd.read_csv(filepath, sep=' ', skiprows=9, names=col_names)
        
    return df



# Function to count number of contacts between the two molekular chains
def CountContacts(type_1_coords, type_2_coords, dist):

    # Compute pairwise distances between type 1 and type 2 particles
    distances = np.zeros((len(type_1_coords), len(type_2_coords)))
    
    for i in range(len(type_1_coords)):
        distances[i,:] = np.sqrt((type_1_coords[i,0]-type_2_coords[:,0])**2 + (type_1_coords[i,1]-type_2_coords[:,1])**2)

    
    count_contacts = np.sum(distances <= dist)
    
    return count_contacts, distances


# Gives the pressure in SI-units (Pa)
def ContactPotential(distances, r_1, acc, E=100*(1e6)):
    # Define constants
    r_2 = 1.5                       # Radius top chain, micrometers
    PR = 0.3                        # Poissons ratio 

    # Calculate effective Youngs modulus
    effective_youngs = ((1-PR**2)/E + (1-PR**2)/E)**(-1)

    # Calculate particle overlap, 0 where there is no overlap
    particle_overlap = np.where(distances <= r_1+r_2, ((r_1+r_2)-distances), 0)

    contact_rows, contact_cols = np.nonzero(particle_overlap) 
   # particle_overlap = np.array([[round(y, 12) for y in x]for x in particle_overlap])

    
    effective_radius = (r_1*r_2)/(r_1+r_2)
    applied_force = 4/3 * effective_youngs * effective_radius**(1/2) * particle_overlap[contact_rows, contact_cols]**(3/2)
    
    # Calculate the maximum contact pressure (pressure in the centre of the contact area)
    maximum_contact_presssure = 1/np.pi *((6 * applied_force * effective_youngs**2) / effective_radius**2)**(1/3)
    
    # Calculate potential in contact areas
    potential = (8/15) * effective_youngs*effective_radius**(1/2) * (particle_overlap)**(5/2)

    return potential, particle_overlap, maximum_contact_presssure, applied_force#, part_energy#,maximum_contact_presssure, contact_radius, particle_overlap[contact_rows, contact_cols], energy, applied_force
 


# Returns 2D array showing placement of contacts
def AllContacts(type_1_coords, type_2_coords, dist):
    conts = np.zeros(len(type_2_coords))
    for i in range(len(type_2_coords)):
        if (np.sqrt((type_1_coords[:,0]-type_2_coords[i,0])**2 + (type_1_coords[:,1]-type_2_coords[i,1])**2) <= dist).any():
            conts[i] = 1
    return conts




# Find potential in system given its dump-file
def find_potential(filename, acc=20, youngs=100*(1e6)):
    #Read data from file
    df = ReadDumpFile(filename) 

    # Separate type 1 and type 2 particles
    type_1_df = df[df.type == 1]
    type_2_df = df[df.type == 2]

    # Extract the relevant coordinates (xu, zu) into numpy arrays
    type_1_coords = type_1_df[['xu', 'zu']].to_numpy()
    type_2_coords = type_2_df[['xu', 'zu']].to_numpy()

    delta_h = np.abs(type_1_coords[0][1] - type_2_coords[0][1])

    r_1 = (df['diameter'][0]) / 2
    dist = 1.5 + r_1
  
    # Count contacts and find distances between all points
    count_contacts, distances = CountContacts(type_1_coords, type_2_coords, dist)
    all_contacts = AllContacts(type_1_coords, type_2_coords, dist)
    
    # Find potential, overlap and maximum contact pressure
    potential, particle_overlap, maximum_contact_pres, applied_force = ContactPotential(distances, r_1, acc, E=youngs)
    #particle_overlap_arr = np.array([round(y, acc) for x in particle_overlap for y in x])
   
    #ocurr_dict = {}
    #for i in range(len(set(particle_overlap_arr))):
        #ocurr_dict[list(set(particle_overlap_arr))[i]] = list(particle_overlap_arr).count(list(set(particle_overlap_arr))[i])
        #print(list(overlaps_arr).count(list(set(overlaps_arr))[i]))

    # Change units to SI-units (from pg/micro s^2) to Pa
    potential *= 10**(-15) #Nm
    maximum_contact_pres *= 10**(3)
    

    return potential, count_contacts, particle_overlap, maximum_contact_pres, all_contacts, delta_h


# Find the distances in the contacts given a dump file
def find_contact_distance(filename, youngs=100*(1e6)):
    #Read data from file
    df = ReadDumpFile(filename) 
    
    # Separate type 1 and type 2 particles
    type_1_df = df[df.type == 1]
    type_2_df = df[df.type == 2]

    # Extract the relevant coordinates (xu, zu) into numpy arrays
    type_1_coords = type_1_df[['xu', 'zu']].to_numpy()
    type_2_coords = type_2_df[['xu', 'zu']].to_numpy()
    
    r_1 = (df['diameter'][0]) / 2
    dist = 1.5 + r_1
  
    # Count contacts and find distances between all points
    count_contacts, distances = CountContacts(type_1_coords, type_2_coords, dist)

    mask = distances <= 3.0
    # Apply the transformations directly on the masked values
    distances_ = np.where(mask, distances, 0)

    return distances_, type_1_coords, type_2_coords



# Find potential when there is only one particletype
def find_potential_one_type(filename, acc=15, a=1,youngs=100*(1e6)):
    #Read data from file
    df = ReadDumpFile(filename) 
    
    # Separate type 1 and type 2 particles
    type_1_df = df[df.id <= a]
    type_2_df = df[df.id > a]

    # Extract the relevant coordinates (xu, zu) into numpy arrays
    type_1_coords = type_1_df[['xu', 'zu']].to_numpy()
    type_2_coords = type_2_df[['xu', 'zu']].to_numpy()

    r_1 = (df['diameter'][0]) / 2
    dist = 1.5 + r_1
  
    # Count contacts and find distances between all points
    count_contacts, distances = CountContacts(type_1_coords, type_2_coords, dist)
    

    # Find maximum contact pressure
    potential, particle_overlap, maximum_contact_presssure, applied_force = ContactPotential(distances, r_1, acc, E=youngs)

    # Change units to SI-units (from pg/micro s^2) to Pa
    potential *= 10**(-15)

    return potential, count_contacts, particle_overlap



# Find potential for case of one type of particle and 3 particles in total
def find_potential_one_type_3(filename, acc=7, youngs=100*(1e6)):
    #Read data from file
    df = ReadDumpFile(filename) 

    # Extract the relevant coordinates (xu, zu) into numpy arrays
    coords = df[['xu', 'zu']].to_numpy()
    type_1_coords = np.array(coords[0:2,:])
    type_2_coords = np.array([coords[2,:]])

    #find_angles(type_1_coords[0], type_2_coords[0])
    #find_angles(type_1_coords[1], type_2_coords[0])
    #find_unit_vectors(type_1_coords[0], type_2_coords[0])
    #find_unit_vectors(type_1_coords[1], type_2_coords[0])
    r_1 = (df['diameter'][0]) / 2
    dist = 1.5 + r_1
  
    # Count contacts and find distances between all points
    count_contacts, distances = CountContacts(type_1_coords, type_2_coords, dist)
    

    # Find maximum contact pressure
    potential, particle_overlap = ContactPotential(distances, r_1, acc, E=youngs)

    # Change units to SI-units (from pg/micro s^2) to Pa
    potential *= 10**(-15)

    return potential, count_contacts, particle_overlap





# Find unit vector beteen particles from coordinates
def find_unit_vectors(type_1_coords_, type_2_coords_):
    a = np.array([type_1_coords_[0]-type_2_coords_[0], type_1_coords_[1]-type_2_coords_[1]])
    a /= np.linalg.norm(a)
    print(a)


############################### NYE FUNKSJONER FRA d_prob #########################################

# Find index for specific value in array 
def find_index(value, array):
    array = np.asarray(array)
    idx = np.argmin((np.abs(array - value)))
    return idx


# Function for circle
def particle(x, r, h=0, k=0):
    return np.sqrt(r**2-(x-h)**2) + k


# Function to fit particle chain
def f_func(x, r, l, amp, h=0):
    return amp * np.sin((2*np.pi*x)/l + np.pi/2)


# Find tangent 
def find_tangent(x1, y1, h=0, k=0):
    dx = y1 - k
    dy = -(x1 - h)

    # Normalize the direction vector to get the unit vector
    magnitude = np.sqrt(dx**2 + dy**2)
    unit_vector = [dx / magnitude, dy / magnitude]
    return unit_vector


# Find the perpendicular to a vector
def perpendicular(a) :
    if len(a) == 0:
        b = [1, 0]
    else:
        b = np.empty_like(a)
        b[0] = -a[1]
        b[1] = a[0]
    return b


# Normalize a vector
def normalize(a):
    a = np.array(a)
    return a/np.linalg.norm(a)


# Find the angle between two vectors
def find_angle(u,v):
    cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return angle_rad


# Find the nearest particle 
def find_nearest(array,value):
    ind = 0 
    min_diff = np.abs(array[0] - value)
 
    for i, elem in enumerate(array):
        diff = np.abs(elem - value)
        if diff < min_diff:
            min_diff = diff
            ind = i

    return ind




## Finding the particle overlap with the three methods

def find_d_from_curve(d_curve, x_t, x_coords_i, x_coords_j, z_i, z_j, effective_radius, acc):
    d_array_og = np.zeros(len(x_t))
    d_array_first = np.zeros(len(x_t))
    d_array_second = np.zeros(len(x_t))
    
    d_temp_array = np.zeros(len(x_t))

    for i in range(len(x_t)):
        # Contact if the function is above zero
        if d_curve[i] >= 0:
            d_temp_array[i] = d_curve[i]
            
            #Onlu enters for the last point in the contact
            if i >= len(x_t)-1 or d_curve[i+1] <= 0:
                max_d = np.max(d_temp_array)
                index_max_d = np.where(d_temp_array == max_d)[0][0]

                # Finds the first and last index of this contact
                first_index_contact = np.argwhere(d_temp_array>0)[0][0]
                last_index_contact = np.argwhere(d_temp_array>0)[-1][0]


                # Find nearest particles to contact
                nearest_index_x_i = find_nearest(x_coords_i, x_t[index_max_d])
                nearest_index_x_j = find_nearest(x_coords_j, x_t[index_max_d])
                
                # No corrections done to d
                d_array_og[i] = round(max_d, acc)
                
                # vector between the particles
                BC = [x_coords_j[nearest_index_x_j] - x_coords_i[nearest_index_x_i], z_j-z_i]

                # Find L (different at ends of the chain)
                if first_index_contact == index_max_d or last_index_contact == index_max_d:
                    L = 2*(x_t[last_index_contact]-x_t[first_index_contact])
                else: 
                    L = x_t[last_index_contact]-x_t[first_index_contact]

                angle = find_angle([0, max_d], BC)
                a = (L / (2 * np.cos(angle)))
                

                # calculate the approximate d 
                #d = effective_radius-np.sqrt(effective_radius**2-a**2)
                d=a**2/(2*effective_radius)
                

                d_array_second[i] = round(d, acc)
                d_array_first[i] = round(np.abs(max_d)*np.cos(angle),acc)

                # Reset temp array
                d_temp_array = np.zeros(len(x_t))
                
    return d_array_og, d_array_first, d_array_second




# Calculate the analytical function
def find_d_curve(filename, E=100*(1e6), PR=0.3, acc=8):  
    effective_youngs = ((1-PR**2)/E + (1-PR**2)/E)**(-1)

    df = ReadDumpFile(filename)
    df_i = df[df.type == 1]
    df_j = df[df.type == 2]

    # Extract the relevant coordinates (xu, zu) into numpy arrays
    coords_i = df_i[['xu', 'zu', 'diameter']].to_numpy()
    coords_j = df_j[['xu', 'zu', 'diameter']].to_numpy()
    x_coords_i = df_i['xu'].to_numpy()
    x_coords_j = df_j['xu'].to_numpy()

    # Get the radius of the particles and calculate the effective radius
    r_i = coords_i[0][2] / 2
    r_j = coords_j[0][2] / 2
    effective_radius = (r_i*r_j)/(r_i+r_j)                  
    
    # Lattice distance (wavelength) 
    l_i = np.abs(coords_i[0][0] - coords_i[1][0])
    l_j = np.abs(coords_j[0][0] - coords_j[1][0])

    # z_positions for the chains
    z_i = coords_i[0][1]
    z_j = coords_j[0][1]

    # Create x-array
    x_t = np.linspace(np.max([coords_i[0][0],coords_j[0][0]])-np.max([r_i, r_j]),np.min([coords_i[-1][0],coords_j[-1][0]])+np.max([r_i, r_j]), len(coords_j)*15000)
    #x_t = np.linspace(np.max([coords_i[0][0],coords_j[0][0]]),np.min([coords_i[-1][0],coords_j[-1][0]]), len(coords_j)*15000)
    
    # amplitude of the sine-waves
    amp_i = np.abs((r_i**2*l_i**2)/(4*np.pi**2*(r_i**2)**(3/2)))
    amp_j = np.abs((r_j**2*l_j**2)/(4*np.pi**2*(r_j**2)**(3/2)))

    # Sine-waved for each chain
    y_i_t = f_func(x_t, r_i, l_i, amp_i)
    y_j_t = f_func(x_t, r_j, l_j, amp_j)

    # Combined sine-wave 
    sine_curve_combined = (z_i+(r_i-amp_i)+ y_i_t)-(z_j-(r_j-amp_j)-y_j_t)

    # Find d from the sine-curve
    d_array_og, d_array_first, d_array_second = find_d_from_curve(sine_curve_combined, x_t, x_coords_i, x_coords_j, z_i, z_j, effective_radius, acc)

    # Calculate contacts and potential from d
    a = effective_youngs * effective_radius**(1/2)
    contacts_og = np.count_nonzero(d_array_og)
    contacts_first = np.count_nonzero(d_array_first)
    contacts_second = np.count_nonzero(d_array_second)

    #max_overlap_og = np.max(d_array_og)
    #max_overlap_first = np.max(d_array_first)
    #max_overlap_second = np.max(d_array_second)
    
    if contacts_og != contacts_first or contacts_og != contacts_second:
        print("ERROR! The number of contacts differs for the different methods.")

    og_potential = (8/15)*a*np.sum((np.array(np.abs(d_array_og)))**(5/2)) *10**(-15)
    first_potential = (8/15)*a*np.sum((np.array(d_array_first))**(5/2)) *10**(-15)
    second_potential = (8/15)*a*np.sum((np.array(d_array_second))**(5/2)) *10**(-15)


    #print('Number of contacts: ', np.count_nonzero(d_array_og))
    #print('No correction for d analytic potential: ', og_potential)
    #print('d: ', set(d_array_og))
    #print('First method analytic potential: ', first_potential)
    #print('d: ', set(d_array_first))
    #print('Second method analytic potential: ', second_potential)
    #print('d: ', set(d_array_first))
    
    return contacts_og, og_potential, first_potential, second_potential, d_array_og[np.where(d_array_og!=0)], d_array_first[np.where(d_array_first!=0)], d_array_second[np.where(d_array_second!=0)]#max_overlap_og, max_overlap_first, max_overlap_second
    

def mse(y_pred, y_acc):
    return (np.sum((y_pred-y_acc)**2))/len(y_pred)

