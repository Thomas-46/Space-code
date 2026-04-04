# Constantes fondamentales

# Constante gravitationnelle universelle
G = 6.67430 * 10**(-11)

# Vitesse de la lumière dans le vide		
c = 299792458299792458

# Masses
M_Soleil = 1.98847 * (10**30)
M_Terre	= 5.97219 * (10**24)
M_Lune = 7.34767309 * (10**22)

# Constantes gravitationnelles
Mu_Soleil = 1.32712440018 * (10**20)
Mu_Terre = 3.986004418 * (10**14)
Mu_Lune = 4.9048695 * (10**12)

# Constantes Mathématiques
PI = 3.141592653589793
TWO_PI = 2 *PI


# ----------------------------------------------------------
# CONSTANTES PHYSIQUES NORMALISÉES
# ----------------------------------------------------------
# Unités : distance en UA, masse en masse solaire, temps en jours
AU = 1.495978707*(10**11)
Day = 86400
G_norm = G * (AU**(-3)) * M_Soleil * (Day**2)


# ----------------------------------------------------------
# CONDITIONS INITIALES DES PLANÈTES (Positions [UA], Vitesses [UA/jour], Masse [M_sol])
# ----------------------------------------------------------

# Soleil
M_Soleil_norm   = 1.0
X0_Soleil_norm  = -1.6132380e-3
Y0_Soleil_norm  =  2.3674938e-3
Z0_Soleil_norm  = -3.4999903e-5
VX0_Soleil_norm =  6.0453208e-6
VY0_Soleil_norm = -1.8980631e-6
VZ0_Soleil_norm = -1.3060634e-7

# Mercure
M_Mercure_norm   = 1.6600e-7
X0_Mercure_norm  =  3.4942761e-1
Y0_Mercure_norm  =  1.0619088e-2
Z0_Mercure_norm  = -3.1182970e-2
VX0_Mercure_norm = -6.4674481e-3
VY0_Mercure_norm =  2.9372684e-2
VZ0_Mercure_norm =  2.9939160e-3

# Vénus
M_Venus_norm   = 2.4476e-6
X0_Venus_norm  = -5.7685710e-1
Y0_Venus_norm  =  4.2639983e-1
Z0_Venus_norm  =  3.9039256e-2
VX0_Venus_norm = -1.2159929e-2
VY0_Venus_norm = -1.6321487e-2
VZ0_Venus_norm =  4.7839154e-4

# Terre
M_Terre_norm   = 3.0032e-6
X0_Terre_norm  =  6.8900355e-1
Y0_Terre_norm  =  7.0799513e-1
Z0_Terre_norm  = -5.4762805e-5
VX0_Terre_norm = -1.2606830e-2
VY0_Terre_norm =  1.1932472e-2
VZ0_Terre_norm = -5.3326343e-7

# Mars
M_Mars_norm   = 3.2268e-7
X0_Mars_norm  =  4.3534005e-1
Y0_Mars_norm  = -1.3548686
Z0_Mars_norm  = -3.9100527e-2
VX0_Mars_norm =  1.3851438e-2
VY0_Mars_norm =  5.5002850e-3
VZ0_Mars_norm = -2.2479471e-4

# Jupiter
M_Jupiter_norm   = 9.5425e-4
X0_Jupiter_norm  =  1.8147606
Y0_Jupiter_norm  =  4.7059241
Z0_Jupiter_norm  = -6.0234633e-2
VX0_Jupiter_norm = -7.1329432e-3
VY0_Jupiter_norm =  3.0767709e-3
VZ0_Jupiter_norm =  1.4683761e-4

# Saturne
M_Saturne_norm   = 2.8572e-4
X0_Saturne_norm  = -8.2319838
Y0_Saturne_norm  = -5.2651142
Z0_Saturne_norm  =  4.1915711e-1
VX0_Saturne_norm =  2.7027945e-3
VY0_Saturne_norm = -4.7129500e-3
VZ0_Saturne_norm = -2.5483556e-5

# Uranus
M_Uranus_norm   = 4.3643e-5
X0_Uranus_norm  =  1.9916762e+01
Y0_Uranus_norm  =  2.3700466
Z0_Uranus_norm  = -2.4922728e-1
VX0_Uranus_norm = -4.9354675e-4
VY0_Uranus_norm =  3.7222065e-3
VZ0_Uranus_norm =  2.0227673e-5

# Neptune
M_Neptune_norm   = 5.1486e-5
X0_Neptune_norm  =  2.6481558e+01
Y0_Neptune_norm  = -1.4076556e+01
Z0_Neptune_norm  = -3.2040642e-1
VX0_Neptune_norm =  1.4527105e-3
VY0_Neptune_norm =  2.7901455e-3
VZ0_Neptune_norm = -9.1342324e-5

# Pluton
M_Pluton_norm   = 6.6060e-9
X0_Pluton_norm  =  4.9371124
Y0_Pluton_norm  = -3.1890424e+01
Z0_Pluton_norm  =  1.9843679
VX0_Pluton_norm =  3.1711953e-3
VY0_Pluton_norm = -1.4100225e-4
VZ0_Pluton_norm = -9.0807325e-4

# Lune (The Moon)
# --- LUNE (Basée sur tes positions Terre) ---
M_Lune_norm   = 3.6942e-8      

# Position : On décale la Lune sur l'axe X par rapport à la Terre (+ 0.00257 UA)
X0_Lune_norm  = 6.8900355e-1 + 0.00257  
Y0_Lune_norm  = 7.0799513e-1
Z0_Lune_norm  = -5.4762805e-5

# Vitesse : On prend la vitesse de la Terre et on ajoute la vitesse orbitale 
# de la Lune sur l'axe Y pour qu'elle soit perpendiculaire au décalage X.
VX0_Lune_norm = -1.2606830e-2
VY0_Lune_norm = 1.1932472e-2 + 0.00059  
VZ0_Lune_norm = -5.3326343e-7

# --- VISUALISATION ---
# --- ÉLÉMENTS KÉPLÉRIENS COMPLETS (Moyennes approximatives) ---
# [a, e, omega (en degrés)]
KEPLER_ELEMENTS = {
    "Mercure": [0.3871, 0.2056, 77.4],
    "Vénus":   [0.7233, 0.0067, 131.5],
    "Terre":   [1.0000, 0.0167, 102.9],
    "Mars":    [1.5237, 0.0934, 336.0],
    "Jupiter": [5.2034, 0.0484, 14.7],
    "Saturne": [9.5371, 0.0541, 92.4],
    "Uranus":  [19.191, 0.0472, 170.9],
    "Neptune": [30.069, 0.0086, 44.9],
    "Pluton":  [39.482, 0.2488, 224.0]
}
