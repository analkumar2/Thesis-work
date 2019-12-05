#exec(open('runmultiple.py').read())


for i in range(10):
    exec(open('RandomModel.py').read())
    if Model['Error']<30:
        print(Model)
        plt.plot(tvec25pA, Vmvec25pA, label='Model25pA')
        plt.plot(expData['25pA'][0], expData['25pA'][1], label='Exp25pA')
        plt.plot(tvec50pA, Vmvec50pA, label='Model50pA')
        plt.plot(expData['50pA'][0], expData['50pA'][1], label='Exp50pA')
        plt.plot(tvec150pA, Vmvec150pA, label='Model150pA')
        plt.plot(expData['150pA'][0], expData['150pA'][1], label='Exp150pA')
        plt.legend()
        plt.show()
