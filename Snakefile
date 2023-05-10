rule results:
    output:
        "src/data/dat_maxDistance_1000_final.txt"
    script:
        "src/scripts/CV_pop_create.py"
