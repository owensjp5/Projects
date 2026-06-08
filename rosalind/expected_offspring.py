################################
# Calculate Expected Offspring #
################################
def expectedOffspring(population):
    total = 0
    for i in range(len(population)):
        if i < 3:
            total += 2*population[i]
        elif i==4:
            total += population[i]
        elif i == 3:
            total += 1.5*population[i]
    return total