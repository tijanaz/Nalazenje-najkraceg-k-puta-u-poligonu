import math, random
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LineString, MultiLineString, GeometryCollection, MultiPoint


# Klasa predstavlja jedno resenje za koji se cuva njegov genetski kod (points = izabranih k-2 tacaka)
# i vrednost funkcije prilagodjenosti.
class Solution:

    def __init__(self, points, unfitness):
        self.points = points
        self.unfitness = unfitness


def GeneticAlgorithm(k, A, B, polygon_points):

    # parametri
    reproduction_size = 300
    generation_size = 1000
    init_generation_size = 5000
    mutation_rate = 0.03
    max_iteration = 15
    tournament_k = 20
    crossover_p = 0.5
    crossover_p_mix = 0.2

    solution_size = k-2         # A i B su fiksni, trazimo ostale tacke

    p = Polygon(polygon_points);

    # Vrednosti izmedju kojih se nalaze tacke poligona
    (min_x, min_y, max_x, max_y) = p.bounds
    # pravougaonik koji sadrzi poligon
    polygon_rect = [(min_x-1, min_y-1), (max_x+1, min_y-1), (max_x+1, max_y+1), (min_x-1, max_y+1), (min_x-1, min_y-1)]

    # Formira se niz koji sadrzi sve celobrojne tacke koje pripadaju poligonu
    all_points = []
    for i in range (int(min_x), int(max_x)):
        for j in range (int(min_y), int(max_y)):
            if (i, j) in polygon_points:
                all_points.extend([(i, j)] * 10)   # dodajemo i temena poligona, vise puta da bi se cesce birala
            if p.contains(Point(i, j)):            # if Point(i, j).within(p):
                all_points.append((i, j))

    ############################## Pomocne funkcije ##########################################
    # incicijalizuje se populacija
    def initial_population():

        population = []
        for i in range(generation_size):
            solution = []
            i = 0
            for j in range(solution_size):
                i = random.randint(i, len(all_points) - solution_size + j)
                solution.append(all_points[i])
            population.append(Solution(solution, unfitness(solution)))
        # vraca se niz elemenata klase Solution
        return population

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


    # sto je resenje gore to je veci unfitness (br_preseka*200 + razdaljina)
    def unfitness(solution):
        fit = 0.0

        # dodju se tacke A i B na pocetak i kraj
        coord = [A]
        for point in solution:
            coord.append(point)
        coord.append(B)
        # formira se linija
        line = LineString(coord)

        # dodaje se duzina linije
        fit = line.length

        # dodaje se broj izlaska linije van poligona*200
        fit += num_intersection(coord)*200
        return fit

    # Funkcija bira reproduction_size resenja iz populacije
    def selection(solutions):

        selected_solutions = [
            selection_tournament_pick_one(solutions)
            for i in range(reproduction_size)]

        return selected_solutions

    # Funkcija vraca reproduction_size najboljih resenja
    def selection_elit(solutions):
        solutions.sort(key=lambda s: s.unfitness)
        return solutions[:reproduction_size]

    # Bira najbolju od (tournament_k) broja jedinki koje se izvlace iz populacije.
    def selection_tournament_pick_one(solutions):
            the_chosen_ones = []
            top_i = None
            for i in range(tournament_k):

                pick = random.randint(0, len(solutions)-1)
                the_chosen_ones.append(solutions[pick])
                if top_i == None or the_chosen_ones[i].unfitness < the_chosen_ones[top_i].unfitness:
                    top_i = i
            return the_chosen_ones[top_i]

    def create_generation(for_reproduction):
        # Od jedinki dobijenih u okviru 'for_reproduction' generise novu generaciju

        new_generation = []
        # cuva se najbolje resenje iz prethodne generacije
        new_generation.append(min(solutions, key=lambda s: s.unfitness))
        # sve dok se ne popuni generacija
        while len(new_generation) < generation_size:
            # Biramo dva nasumicno i vrsimo

            parents = random.sample(for_reproduction, 2)
            child1, child2 = crossover_uniform(parents[0].points, parents[1].points)

            # Vrsi se mutacija nakon ukrstanja
            child1 = mutation(child1)
            child2 = mutation(child2)

            # Dodajemo nove hromozome u novu generaciju
            new_generation.append(Solution(child1, unfitness(child1)))
            new_generation.append(Solution(child2, unfitness(child2)))

        return new_generation

    def crossover(a, b):
        # Vrsi jednopoziciono ukrstanje po nasumicnoj poziciji unutar hromozoma.
        cross_point = random.randint(0, len(a)-1)
        ab = a[:cross_point] + b[cross_point:]
        ba = b[:cross_point] + a[cross_point:]
        return (ab, ba)

    def crossover_uniform(a, b):
        # Uniformno ukrstanje po verovatnoci self._crossover_p.
        ab = [(0, 0)] * len(a)
        ba = [(0, 0)] * len(a)
        for i in range(len(a)):
            p1 = random.random()
            p2 = random.random()

            if p1 < crossover_p:
                ab[i] = a[i]
                ba[i] = b[i]
                # Ako je vrednost p2 manja od crossover_p mesaju se i x i y koordinate, a ne samo cele tacke
                if p2 < crossover_p:
                    ab[i] = (b[i][0], ab[i][1])
                    ba[i] = (a[i][0], ba[i][1])
            else:
                ab[i] = b[i]
                ba[i] = a[i]
                if p2 < crossover_p:
                    ab[i] = (ab[i][0], a[i][1])
                    ba[i] = (ba[i][0], b[i][1])
                # temena su cesto uokviru najboljeg resenja pa ih malo cesce biramo
                if a[i] in polygon_points:
                    ab[i] = a[i]
                    ba[i] = a[i]
                if b[i] in polygon_points:
                    ab[i] = b[i]
                    ba[i] = b[i]
        ab.sort(key=lambda x: x[0])
        ba.sort(key=lambda x: x[0])
        return (ab, ba)

    def mutation(solution):
        # Vrsi mutaciju nad resenjem (hromozomom) sa verovatnocom mutation_rate
        t = random.random()
        if t < mutation_rate:
            i = random.randint(0, len(solution)-1)
            solution[i] = random.choice(all_points)
        return solution

    ################################### Glavni deo algoritma #################################################

    # Generise se pocetna populaciju jedinki (resenja) i racuna se prilagodjenost svake jedinke u populaciji
    solutions = initial_population()

    # pamti se koliko je poslednjih najboljih resenja imalo istu funkciju prilagodjenosti
    number_last_same = 0
    last_unfit = -1

    # Sve dok uslov zaustavljanja nije zadovoljen
    current_iteration = 0

    while current_iteration != max_iteration:

        # Prikaz trenutnog stanja algoritma
        print("----------------------------------------------------------")
        print("Iteration: %d" % current_iteration)

        the_sum = sum(solution.unfitness for solution in solutions)
        print("Reproduction chromos sum unfitness: %d" % the_sum)

        top_unfit = min(solutions, key=lambda s: s.unfitness).unfitness
        print("Top solution: %f" % top_unfit)
        print("----------------------------------------------------------")

        # Bira se skup jedinki za reprodukciju iz populacije
        for_reproduction = selection(solutions)

        # Primenom operatora ukrstanja i mutacije kreiraju se nove jedinke i racuna se njihova prilagodjenost.
        # Dobijene jedinke predstavljaju novu generaciju.
        solutions = create_generation(for_reproduction)

        if top_unfit == last_unfit:
            number_last_same += 1
        else:
            number_last_same = 0
        # ako u poslednjih 20 iteracija nema promene prekida se algoritam
        if number_last_same > 40:
            # prikazuje se rezultat i vraca se najbolje resenje
            top_solution = min(solutions, key=lambda chromo: chromo.unfitness)
            print("Broj preseka: ")
            print(num_intersection(top_solution.points))
            return top_solution

        # Prelazak u sledecu iteraciju
        current_iteration += 1

    # Ako je dostignut maksimalan broj iteracija prikazuje se rezultat i vraca se najbolje resenje
    top_solution = min(solutions, key=lambda chromo: chromo.unfitness)
    print("Broj preseka: ")
    print(num_intersection(top_solution.points))
    return top_solution

#######################################################################################################

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
    #A = (160, 321)
    #B = (234, 144)
    
    #2.
    #A = (362.55631868,  220.7947158)
    #B = (260.32511082,  151.84525141)
    
    #3.
    A = (221.70344421,  230.73104356)
    B = (159.60986738,  323.70420134)
    
    top_solution = GeneticAlgorithm(5, A, B, polygon_points)
    print("Solution: %s unfitness: %d" % (top_solution.points, top_solution.unfitness))

    coord = [A]
    for point in top_solution.points:
        coord.append(point)
    coord.append(B)


    plt.figure()
    plt.plot([p[0] for p in coord], [p[1] for p in coord], c = 'r')

    polygon_points_x = [p[0] for p in polygon_points]
    polygon_points_y = [p[1] for p in polygon_points]
    polygon_points_x.append(polygon_points_x[0])
    polygon_points_y.append(polygon_points_y[0])

    plt.plot(polygon_points_x, polygon_points_y)
    plt.scatter([A[0], B[0]], [A[1], B[1]], s=10, c='g')

    plt.show()
