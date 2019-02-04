import math, random
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LineString, MultiLineString, GeometryCollection, MultiPoint


def PSO(k, A, B, polygon_points):

    solution_dim = (k-2)*2

    # parametri
    num_particles = 150
    max_iteration = 200
    velocity_par = 70

    err_best_g = -1                   # najbolja greska za grupu
    pos_best_g = []                   # najbolje pozicija za grupu

    # pocetne tacke se biraju na liniji AB u jednakim razmacima
    dist_x = 1.0*(B[0] - A[0])/(k-1)
    dist_y = 1.0*(B[1] - A[1])/(k-1)

    inital_dots = []
    for i in range(k-2):
        inital_dots.append(A[0]+dist_x*(i+1))
        inital_dots.append(A[1]+dist_y*(i+1))

    p = Polygon(polygon_points);

    # Vrednosti izmedju kojih se nalaze tacke poligona
    (min_x, min_y, max_x, max_y) = p.bounds

    # donje i gornje granice
    bounds=[]
    for i in range(k-2):
        bounds.extend([(min_x, max_x),(min_y,max_y)])

    # pravougaonik koji sadrzi poligon
    polygon_rect = [(min_x-1, min_y-1), (max_x+1, min_y-1), (max_x+1, max_y+1), (min_x-1, max_y+1), (min_x-1, min_y-1)]


    ################## Pomocne funkcije #################

    # vraca koliko puta linija izlazi van poligona
    def num_intersection(coord):
        # deo linije koji je van poligona se trazi kao presek poligona koji je okvir oko naseg poligona i linije
        out = Polygon(polygon_rect, [polygon_points]).intersection(LineString(coord))

        if out:     # ako linija izlazi van poligona
            if type(out) is MultiLineString:    # ako izlazi iz poligona na vise mesta
                return len(out)
            elif type(out) is GeometryCollection:
                sum = 0
                for o in out:
                    if type(o) is LineString:
                        sum+= 1
                return sum
            elif type(out) is LineString:
                return 1
            elif type(out) is Point:
                return 0
            elif type(out) is MultiPoint:
                return 0
            else:
                return 0
        else:
            return 0



    # funkcija koju pokusavamo da minimizujemo
    # vraca duzinu puta + dodatne vrednosti ako put nije ceo u poligonu
    def func1(solution):
        total=0

        zipped = zip(solution[0::2], solution[1::2])

        coord = [A]
        for point in zipped:
            # ako je neka tacka van poligona, resenje je lose
            if not Polygon(polygon_points).contains(Point(point[0], point[1])):
                return solution_dim*400
            coord.append(point)
        coord.append(B)
        # formira se linija
        line = LineString(coord)

        # dodaje se duzina linije
        total = line.length

        # dodaje se broj izlaska linije van poligona*200
        total += num_intersection(coord)*200
        return total

    ##############################################################

    # Klasa predstavlja cesticu (resenje)
    class Particle:
        def __init__(self, inital_dots):
            self.position_i = []          # pozicija cestice
            self.velocity_i = []          # brzina cestice
            self.pos_best_i = []          # najbolja pozicija
            self.err_best_i = -1          # najbolja greska
            self.err_i      = -1          # pojedinacna greska

            for i in range(0, solution_dim):
                # pocetne brzine se biraju nasumicno
                self.velocity_i.append(random.uniform(-velocity_par, velocity_par))
                self.position_i.append(inital_dots[i])

        # izracunati trenutnu prilagodjenost
        def evaluate(self, costFunc):
            self.err_i = costFunc(self.position_i)

            # proveriti da li je trenutna pozicija najbolja pozicija
            if self.err_i < self.err_best_i or self.err_best_i == -1:

                self.pos_best_i = self.position_i
                self.err_best_i = self.err_i


        # azuriranje brzine cestice
        def update_velocity(self,pos_best_g):
            w  = 0.5       # koliko utice prethodna brzina
            c1 = 1
            c2 = 2

            for i in range(0,solution_dim):
                r1 = random.random()
                r2 = random.random()

                vel_cognitive = c1 * r1 * (self.pos_best_i[i] - self.position_i[i])
                vel_social    = c2 * r2 * (pos_best_g[i] - self.position_i[i])

                self.velocity_i[i] = w * self.velocity_i[i] + vel_cognitive + vel_social

        # azuriranje pozicije cestice na osnovu azurirane brzine
        def update_position(self, bounds):
            for i in range(0, solution_dim):
                self.position_i[i]=self.position_i[i]+self.velocity_i[i]

                # ako je neophodno prilagoditi do maksimalne pozicije
                if self.position_i[i] > bounds[i][1]:
                    self.position_i[i] = bounds[i][1]

                # ako je neophodno prilagoditi do minimalne pozicije
                if self.position_i[i] < bounds[i][0]:
                    self.position_i[i]=bounds[i][0]

    #######################################################


    # ako se linija AB nalazi unutar poligona vraca se AB (od k tacaka)
    if num_intersection([A, B]) == 0:
        return inital_dots

    # uspostavlja se roj
    swarm = []
    for i in range(0, num_particles):
        swarm.append(Particle(inital_dots))

    # glavni deo algoritma
    i=0
    # promenjive za pracenje u koliko poslednjih iteracija nije bilo promene
    number_last_same = 0
    last_err_best_g = -1

    while i < max_iteration:

        # prikaz stanja
        print i
        print pos_best_g
        print err_best_g
        print "-------------------------------------------------"

        # obilazak cestica u roju i azuriranje prilagodjenosti
        for j in range(0, num_particles):
            swarm[j].evaluate(func1)

            # utvrdjuje se da li je trenutna cestica i najbolja
            if swarm[j].err_i < err_best_g or err_best_g == -1:
                pos_best_g = list (swarm[j].position_i)
                err_best_g = float(swarm[j].err_i)

        # obilazak cestica u roju i azuriranje brzina i pozicija
        for j in range(0,num_particles):
            swarm[j].update_velocity(pos_best_g)
            swarm[j].update_position(bounds)

        # azurira se broj poslednjih iteracija u kojima nije doslo do promene
        if err_best_g == last_err_best_g:
            number_last_same += 1
        else:
            number_last_same = 0
            last_err_best_g = err_best_g

        # ako u poslednjih 20 iteracija nema promene prekida se algoritam
        if number_last_same > 20:
            # prikazuje se rezultat i vraca se najbolje resenje
            print 'Rezultat:'
            print pos_best_g
            print err_best_g
            return pos_best_g

        # prelazak u narednu iteraciju
        i+=1

    # prikazuje se rezultat i vraca se najbolje resenje
    print 'Rezultat:'
    print pos_best_g
    print err_best_g
    return pos_best_g


if __name__ == "__main__":

    random.seed()

    #zadat poligon
    polygon_points = [(105, 194), (162, 242), (153, 264), (171, 276), (153, 310),
          (142, 360), (193, 328),
          (250, 248), (233, 326), (256, 343), (288, 378), (277, 307),
          (331, 362), (329, 320), (339, 297), (328, 276), (340, 266),
          (316, 245), (382, 226), (443, 222), (314, 210), (300, 195),
          (267, 225), (319, 125), (278, 160), (223, 130), (262, 192),
          (210,220), (150, 201), (143,197)]
    #1.
    A = (160, 321)
    B = (234, 144)

    #2.
    #A = (362.55631868,  220.7947158)
    #B = (260.32511082,  151.84525141)

    #3.
    #A = (221.70344421,  230.73104356)
    #B = (159.60986738,  323.70420134)

    res = PSO(6, A, B, polygon_points)

    res = zip(res[0::2], res[1::2])

    coord = [A]
    for point in res:
        coord.append(point)
    coord.append(B)

    print(coord)

    plt.figure()
    plt.plot([p[0] for p in coord], [p[1] for p in coord], c = 'r')

    polygon_points_x = [p[0] for p in polygon_points]
    polygon_points_y = [p[1] for p in polygon_points]
    polygon_points_x.append(polygon_points_x[0])
    polygon_points_y.append(polygon_points_y[0])

    plt.plot(polygon_points_x, polygon_points_y)
    plt.scatter([A[0], B[0]], [A[1], B[1]], s=10, c='g')

    plt.show()
