from exercise2 import A2_1_Kainz_Jakob, A2_2_Kainz_Jakob
from exercise3 import A3_Kainz_Jakob_intern
from exercise4 import A4_Kainz_Jakob_intern, config


isValid = False
while not isValid:
    ex_number = input('Which task shall be done? Use numbers to select exercise (Currently available tasks: 2, 3, 4) \n')
    if ex_number.isdigit():
        if ex_number == "2":
            ex_number = input('2nd exercise has two tasks, which task shall be done? Use numerics 1 for dna to protein, 2 for GenBank exercise')
            if ex_number.isdigit():
                if ex_number == "1":
                    isValid = True
                    A2_1_Kainz_Jakob.translate_dna_to_protein()
                elif ex_number == "2":
                    isValid = True
                    A2_2_Kainz_Jakob.write_records()
        elif ex_number == "3":
            isValid = True
            A3_Kainz_Jakob_intern.kmeans_initial()
        elif ex_number == "4":
            isValid = True
            A4_Kainz_Jakob_intern.create_tree(config.ENTRAZ_EMAIL, config.GENE_NAME_ALTERNATIVE_3)

    else:
        print('Congrats, you failed to give a valid input, try again')
