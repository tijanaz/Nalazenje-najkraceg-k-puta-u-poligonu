import math, random
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LineString

#kod je preuzet sa interneta
def generatePolygon( ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts ):
    irregularity = clip( irregularity, 0,1 ) * 2*math.pi / numVerts
    spikeyness = clip( spikeyness, 0,1 ) * aveRadius

    # generate n angle steps
    angleSteps = []
    lower = (2*math.pi / numVerts) - irregularity
    upper = (2*math.pi / numVerts) + irregularity
    sum = 0
    for i in range(numVerts) :
        tmp = random.uniform(lower, upper)
        angleSteps.append( tmp )
        sum = sum + tmp

    # normalize the steps so that point 0 and point n+1 are the same
    k = sum / (2*math.pi)
    for i in range(numVerts) :
        angleSteps[i] = angleSteps[i] / k

    # now generate the points
    points = []
    angle = random.uniform(0, 2*math.pi)
    for i in range(numVerts) :
        r_i = clip( random.gauss(aveRadius, spikeyness), 0, 2*aveRadius )
        x = ctrX + r_i*math.cos(angle)
        y = ctrY + r_i*math.sin(angle)
        points.append( (int(x),int(y)) )

        angle = angle + angleSteps[i]

    return points

def clip(x, min, max):
    if( min > max ) :  return x
    elif( x < min ) :  return min
    elif( x > max ) :  return max
    else :             return x

def random_point_within(poly):
    min_x, min_y, max_x, max_y = poly.bounds

    x = random.uniform(min_x, max_x)
    x_line = LineString([(x, min_y), (x, max_y)])
    x_line_intercept_min, x_line_intercept_max = x_line.intersection(poly).xy[1].tolist()
    y = random.uniform(x_line_intercept_min, x_line_intercept_max)

    return Point([x, y])

def get_random_point_in_polygon(poly):
     (minx, miny, maxx, maxy) = poly.bounds
     while True:
         p = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
         if poly.contains(p):
             return p
         
#vraca broj presecnih tacaka poliona i linije
#ako presecnih tacaka ima vise, vraca i koord presecne tacke i redni broj ivice
#poligona gde je presek
def intersect(x, y, line, a, b, second):
    intersect_count = 0
    intersect_line = 0
    ind = 0
    ind_counter = 1
    line2 = line;
    for i in range (0, 30):
        ind_counter = 1
        
        if i != 29:
            line1 = LineString([(x[i], y[i]), (x[i+1], y[i+1])])
        else:
            line1 = LineString([(x[i], y[i]), (x[0], y[0])])
    
        p = line1.intersection(line2)
        if p:
            #print("true")
                    
            for j in range (0, 30):
                if (p.x == x[j] and p.y == y[j]):
                    ind_counter = 0
            if ind_counter==1:
                intersect_count = intersect_count + 1
            if ind == 0 and ind_counter == 1:
                p_x = p.x
                p_y = p.y
                ind = 1
                intersect_line = i
            elif  second == 1 and ind == 1 and ind_counter == 1:
                p_x = p.x
                p_y = p.y
                ind = 2
                intersect_line = i
            
                
                
    if intersect_count >=1:
        return [intersect_count, p_x, p_y, intersect_line]
    else:
        return [intersect_count]

#a_x i a-y - koord presecne tacke 
#centar poligona podesen na (250, 250)
#vraca koordinate najblizeg temena od presecne tacke, koja je najbliza centru
def nearest_vertice(n, m, a_x, a_y, intersect_line):
    
    #rastojanje od levog temena (u odnosu na presecnu) do centra
    left_to_center = math.sqrt((n[intersect_line]-250)*(n[intersect_line]-250) + (m[intersect_line]-250)*(m[intersect_line]-250))
    
    #rastojanje od desnog temena (u odnosu na presecnu) do centra
    right_to_center = math.sqrt((n[intersect_line+1]-250)*(n[intersect_line+1]-250) + (m[intersect_line+1]-250)*(m[intersect_line+1]-250))
    
    if left_to_center < right_to_center:
        return [n[intersect_line], m[intersect_line]]
    else:
         return [n[intersect_line+1], m[intersect_line+1]]
     
def algorytm(x, y, lineAB,a ,b, n, m):
    
    #pravljenje niza (1 za x koord, 1 za y koord) gde ce se cuvati temena konacnog puta
    #tipa - > list
    array_vertices_x = []
    array_vertices_y = []
    array_vertices_x.append(a[0])
    array_vertices_x.append(b[0])
    array_vertices_y.append(a[1])
    array_vertices_y.append(b[1])
    
    second = 0
    k = 1;
    counter = intersect(x,y,lineAB, a, b, second)[0]
    length = 0
    pointer = 1
    while counter >=1 :
        
        counter = intersect(x,y,lineAB, a, b, second)[0]
        
        #print("----Koordinate presecne tacke:")
        #print(intersect(x,y,lineAB, a, b, second)[1])
        #print(intersect(x,y,lineAB, a, b, second)[2])
        
        #nova pocetna tacka za nastvak puta, nadjeno najblize teme
        p = nearest_vertice(n, m, intersect(x,y,lineAB, a, b, second)[1], intersect(x,y,lineAB, a, b, second)[2], intersect(x,y,lineAB, a, b, second)[3])
        #print("----Najblize teme:")
        #print(p)
        
        #provera da li je nova pocetna tacka bliza pocetnoj ili krajnjoj
        #ako je bliza pocetnoj-nadjeno teme je nova pocetna
        #a ako je bliza krajnjoj-onda je nadjeo teme nova krajnja
        a_to_p = math.sqrt((p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]))
        b_to_p = math.sqrt((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]))
        
        #menjanje pocetne ili krajnje tacke
        if a_to_p < b_to_p:
            #print("----Novi pocetak:")
            plt.plot([p[0], a[0]], [p[1], a[1]])
            length = length + math.sqrt((p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]))
            a = p
            #print(a)
            if a[0] > b[0]:
                pom = a
                a = b
                b = pom
            if len(array_vertices_x) == 2:
                array_vertices_x.insert(pointer, p[0])
                array_vertices_y.insert(pointer, p[1])
                pointer = pointer + 1
            else:
                array_vertices_x.insert(pointer, p[0])
                array_vertices_y.insert(pointer, p[1])
                pointer = pointer + 1
        else:
            #print("----Novi kraj:")
            plt.plot([p[0], b[0]], [p[1], b[1]])
            length = length + math.sqrt((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]))
            b = p
            #print(b)
            if a[0] > b[0]:
                pom = a
                a = b
                b = pom
            array_vertices_x.insert(pointer, p[0])
            array_vertices_y.insert(pointer, p[1])
        
        
        
        #nova linija
        lineAB = LineString([(a[0], a[1]), (b[0], b[1])])
        second = 1
        counter = intersect(x,y,lineAB, a, b, second)[0]
        
        plt.scatter(a[0], a[1],s=10, c='g')
        plt.scatter(b[0], b[1],s=10, c='g')
        plt.plot([a[0], b[0]], [a[1], b[1]], 'b')
        
        k = k+1
        
    #print("----Broj duzi(k) pre:")
    #print(k)

    length = length + math.sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
    #print("----Duzina putanje pre:")
    #print(length) 

    return [k, length, array_vertices_x, array_vertices_y]

#trazi da li je moguca putanja izmedju 2 nesusdne tacke cija je duzina manja
def algorytm2(length, array_vertices_x, array_vertices_y, a, b, x, y, k):
    length_new = 0
    if k>2:
        for i in range (0, len(array_vertices_x)-2):
            if k > 2:
                a = array_vertices_x
                b = array_vertices_y
                length_new = math.sqrt((a[i]-b[i])*(a[i]-b[i]) + (a[i+2]-b[i+2])*(a[i+2]-b[i+2]))
                #print("Nova duzina:")
                #print(length_new)
                line = LineString([(a[i], b[i]), (a[i+2], b[i+2])])
                length_two_lines = math.sqrt((a[i]-b[i])*(a[i]-b[i]) + (a[i+1]-b[i+1])*(a[i+1]-b[i+1])) + math.sqrt((a[i+1]-b[i+1])*(a[i+1]-b[i+1]) + (a[i+2]-b[i+2])*(a[i+2]-b[i+2]))
                t = array_vertices_x[i+1]
                #print("stara duzina:")
                #print(length_two_lines)
                #ako nema preseka i nova duzina je manja
                if (intersect(x,y,line, a, b, 1)[0] == 0) and (length_new < length_two_lines):
                    plt.plot([a[i], a[i+2]], [b[i], b[i+2]])
                    k= k-1
                    length = length-length_two_lines
                    length = length + length_new
                    #azuriramo listu temena putanje
                    array_x_tmp = []
                    array_y_tmp = []
                    for j in range (0, len(array_vertices_x)):
                        if array_vertices_x[j] != t:
                            array_x_tmp.append(array_vertices_x[j])
                            array_y_tmp.append(array_vertices_y[j])
                    #print("Nova temena:")
                    #print(array_x_tmp)
                    array_vertices_x = array_x_tmp
                    #print("Nova temena:")
                    #print(array_vertices_x)
                    array_vertices_y = array_y_tmp
                    
    #print("*Nova temena:")
    #print(array_vertices_x)
    return [k, length, array_vertices_x, array_vertices_y]

def k_algorytm(max_k, k, array_vertices_x, array_vertices_y, x, y, length):
    
    list_new_array_vertices_x = []
    list_new_array_vertices_y = []
    list_new_length = []
    
    for i in range(0, len(array_vertices_x)-2):
        line1 = LineString([(array_vertices_x[i], array_vertices_y[i]), (array_vertices_x[i+1], array_vertices_y[i+1])])
        coef1 = (array_vertices_y[i+1]-array_vertices_y[i])/(array_vertices_x[i+1]-array_vertices_x[i])
        n1 = array_vertices_y[i]-coef1*array_vertices_x[i]
        
        for j in range(i+1, len(array_vertices_x)-1):
            line2 = LineString([(array_vertices_x[j], array_vertices_y[j]), (array_vertices_x[j+1], array_vertices_y[j+1])])
            coef2 = (array_vertices_y[j+1]-array_vertices_y[j])/(array_vertices_x[j+1]-array_vertices_x[j])
            n2 = array_vertices_y[j]-coef2*array_vertices_x[j]
            xy = np.linalg.solve(np.array([[-coef1, 1], [-coef2, 1]]), np.array([n1, n2]))
            #print("Nadjeni presek")
            #print(xy)

            indicator = 1
            #provera da li je presecna tacka xy teme poligona
            for l in range(0, 30):
                if xy.item(0) == x[l] and xy.item(1) == y[l]:
                    print("Presek je teme poligona")
                    indicator = 0
            
            ##TODO provera da li tacka xy pripada poligonu
            if 1:
                if indicator:
                    new_vertices_x = []
                    new_vertices_y = []
                    new_vertices_x = array_vertices_x
                    new_vertices_y = array_vertices_y
                    #print(new_vertices_x)
                    za_izbacvanje_x = []
                    za_izbacvanje_y = []
                    za_izbacvanje_x = new_vertices_x[i+1:j+1]
                    za_izbacvanje_y = new_vertices_y[i+1:j+1]
                    print("Za izbacivanje")
                    print(za_izbacvanje_x)
                    for h in range(0, len(za_izbacvanje_x)):
                        new_vertices_x.remove(za_izbacvanje_x[h])
                        new_vertices_y.remove(za_izbacvanje_y[h])
                    
                    
                    ##print(type(new_vertices_x))
                    ##print(new_vertices_x)
                    new_vertices_x.insert(i+1, xy.item(0))
                    new_vertices_y.insert(i+1, xy.item(1))
                    
                    
                    print(new_vertices_x)
                    
                    #duzina nove moguce putanje
                    new_length = 0
                    for p in range(0, len(new_vertices_x)-1): 
                        new_length = new_length + math.sqrt((new_vertices_x[p]-new_vertices_x[p+1])*(new_vertices_x[p]-new_vertices_x[p+1]) + (new_vertices_y[p]-new_vertices_y[p+1])*(new_vertices_y[p]-new_vertices_y[p+1]))
                        
                        if new_length < length:
                            list_new_array_vertices_x.append(new_vertices_x)
                            list_new_array_vertices_y.append(new_vertices_y)
                            list_new_length.append(new_length)
                            
                        elif max_k < k:
                            list_new_array_vertices_x.append(new_vertices_x)
                            list_new_array_vertices_y.append(new_vertices_y)
                            list_new_length.append(new_length)
                            
    min_len = min(list_new_length)
    if min_len < length:
        index =  list_new_length.index(min_len)
        
        array_vertices_x = list_new_array_vertices_x[index]
        array_vertices_x = list_new_array_vertices_x[index]
        k = k-1
        
    return [array_vertices_x, array_vertices_y, k]

#points = generatePolygon( ctrX=250, ctrY=250, aveRadius=100, irregularity=0.5, spikeyness=0.4, numVerts=30 )

#points = [(97, 198), (154, 230), (125, 240), (111, 265), (152, 239),
          #(120, 250), (183, 275),
          #(175, 309), (209, 308), (205, 341), (235, 379), (258, 372),
          #(255, 270), (290, 307), (325, 327), (318, 296), (355, 289),
          #(371, 273), (339, 242), (420, 202), (344, 207), (278, 229),
          #(321, 171), (329, 123), (267, 170), (245, 155), (236, 192),
          #(190,126), (179, 179), (148,195)]
          
points = [(105, 194), (162, 242), (153, 264), (171, 276), (153, 310),
          (142, 360), (193, 328),
          (250, 248), (233, 326), (256, 343), (288, 378), (277, 307),
          (331, 362), (329, 320), (339, 297), (328, 276), (340, 266),
          (316, 245), (382, 226), (443, 222), (314, 210), (300, 195),
          (267, 225), (319, 125), (278, 160), (223, 130), (262, 192),
          (210,220), (150, 201), (143,197)]

#koordinate temena poligona
x = [p[0] for p in points]
y = [p[1] for p in points]
x.append(x[0])
y.append(y[0])

plt.figure()
plt.plot(x, y)

p = Polygon(points)

#point_A = (179, 327)
#point_B = (325, 343)

#1.
#point_A = Point(160, 321)
#point_B = Point(234, 144)

#2.
#point_A = Point(362.55631868,  220.7947158)
#point_B = Point(260.32511082,  151.84525141)

#3.
point_A = Point(221.70344421,  230.73104356)
point_B = Point(159.60986738,  323.70420134)

n, m =p.exterior.coords.xy #n, m-koordinate temena poligona(tipa array)
#point_A = get_random_point_in_polygon(p)
print(type(point_A))
a = np.array(point_A) #a-koord pocetne tacke (tipa array)
#point_B = get_random_point_in_polygon(p)
b = np.array(point_B) #b-koord krajnje tacke (tipa array)



print("Pocetna tacka:")
print(a)
print("Krajnja tacka:")
print(b)

#zamena pocetne i krajnje tacke, ako treba
if a[0] > b[0]:
    pom = a
    a = b
    b = pom

#pocetna linija, od pocetne do krajnje tacke
lineAB = LineString([(point_A.x, point_A.y), (point_B.x, point_B.y)])
plt.plot([point_A.x, point_B.x], [point_A.y, point_B.y])

max_k = 2

array_vertices_x  = []
array_vertices_y  = []


k = algorytm(x ,y ,lineAB, a, b, n, m)[0]
length = algorytm(x ,y ,lineAB, a, b, n, m)[1]
array_vertices_x = algorytm(x ,y ,lineAB, a, b, n, m)[2]
array_vertices_y = algorytm(x ,y ,lineAB, a, b, n, m)[3]

#print("temena 1 :")
#print(array_vertices_x)
#print(array_vertices_y)


array_vertices_x1 = array_vertices_x
array_vertices_y1 = array_vertices_y
array_vertices_x = algorytm2(length, array_vertices_x1, array_vertices_y1, a, b, x, y, k)[2]
array_vertices_y = algorytm2(length, array_vertices_x1, array_vertices_y1, a, b, x, y, k)[3]
length1 = length
length= algorytm2(length, array_vertices_x1, array_vertices_y1, a, b, x, y, k)[1]
k = algorytm2(length1, array_vertices_x1, array_vertices_y1, a, b, x, y, k)[0]
#print(array_vertices_x)
#print(array_vertices_y)

if max_k < k:
    array_vertices_x1 = array_vertices_x
    array_vertices_y1 = array_vertices_y
    array_vertices_x = k_algorytm(max_k, k, array_vertices_x1, array_vertices_y1, x, y, length)[0]
    array_vertices_y = k_algorytm(max_k, k, array_vertices_x1, array_vertices_y1, x, y, length)[1]
    k = k_algorytm(max_k, k, array_vertices_x1, array_vertices_y1, x, y, length)[2]
    
#k = algorytm2(length, array_vertices_x, array_vertices_y, a, b, x, y, k)[0]
#length = algorytm2(length, array_vertices_x, array_vertices_y, a, b, x, y, k)[1]
#array_vertices_x = algorytm2(length, array_vertices_x, array_vertices_y, a, b, x, y, k)[2]
#array_vertices_y = algorytm2(length, array_vertices_x, array_vertices_y, a, b, x, y, k)[3]

print("Broj k duzi:")
print(k)
print("Duzina linije:")
print(length)
print("Konacna temena:")
print(array_vertices_x)
print(array_vertices_y)

#plt.figure()
plt.plot([p for p in array_vertices_x], [p for p in array_vertices_y], c = 'r')

if max_k < k:
    print("Ne moze se iscrtati putanja uz pomoc max_k duzi")

plt.scatter(point_A.x, point_A.y,s=10, c='m')
plt.scatter(point_B.x, point_B.y,s=10, c='m')
plt.show()