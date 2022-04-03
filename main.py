from exercise2 import A2_1_Kainz_Jakob, A2_2_Kainz_Jakob

isValid = False
while not isValid:
    ex_number = input('Which task shall be done? Use numerics 1 for dna to protein, 2 for GenBank exercise')
    if ex_number.isdigit():
        if ex_number == "1":
            isValid = True
            A2_1_Kainz_Jakob.translate_dna_to_protein()
        elif ex_number == "2":
            isValid = True
            A2_2_Kainz_Jakob.write_records()
    else:
        print('Congrats, you failed to give a valid input, try again')
