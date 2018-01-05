def periodic_p(current,first,L): #for periodic with reference to head bead
    if current - first > 0.5*L:
        return current - L
    elif current - first < -0.5*L:
        return current + L
    else:
        return current

def periodic(dr,L): #give the displacement and box size. Both floats
    if dr > 0.5*L:
        return dr - L
    elif dr < -0.5*L:
        return dr + L
    else:
        return dr