def calc_prop_at_droplet_vapor_layer(prop_droplet, prop_gas, alpha=1/3):
    prop_mix = prop_droplet + alpha*(prop_gas - prop_droplet)
    return prop_mix

def calc_fRe(Re):
    fRe = max(1, Re**0.077)
    return fRe
