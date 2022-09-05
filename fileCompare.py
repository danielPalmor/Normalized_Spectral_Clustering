import numpy as np
import pandas as pd


def main():
    goals = ["wam", "ddg", "lnorm", "jacobi", "spk"]
    k = [0, 3, 0, 2, 0, 3, 0, 4, 0, 10, 0, 5, 0, 20, 0, 7, 0, 10, 0, 19, 0, 25, 0, 30, 0, 10, 40, 0, 2, 0, 2, 0, 5]
    c_or_py = "c"

    twelve = 12
    k_index = 0
    input_num = 0
    c_finished = False
    index = 0
    mod = 1
    num_of_different_files = 0
    output = open("/Users/danielpalmor/Desktop/final_project/differences.txt", "w")
    while index < len(goals):
        try:
            if goals[index] == "jacobi":
                path1 = "/Users/danielpalmor/CLionProjects/final_project/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + ".txt"
                if c_finished:
                    path2 = "/Users/danielpalmor/CLionProjects/final_project/tests_inputs_outputs/outputs/" \
                            + goals[index] + "/" + goals[index] + "_python_sym_matrix_input_" + str(input_num) + ".txt"
                else:
                    path2 = "/Users/danielpalmor/CLionProjects/final_project/tests_inputs_outputs/outputs/" \
                            + goals[index] + "/" + goals[index] + "_" + c_or_py + "_sym_matrix_input_" + str(input_num) \
                            + ".txt"

            elif goals[index] != "spk":
                path1 = "/Users/danielpalmor/CLionProjects/final_project/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + ".txt"

                if c_finished:
                    path2 = "/Users/danielpalmor/CLionProjects/final_project/tests_inputs_outputs/outputs/" \
                            + goals[index] + "/" + goals[index] + "_python_input_" + str(input_num) + ".txt"
                else:
                    path2 = "/Users/danielpalmor/CLionProjects/final_project/tests_inputs_outputs/outputs/" \
                            + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + ".txt"
            else:
                c_finished = True
                c_or_py = "py"
                if k_index >= len(k):
                    return FileNotFoundError
                path1 = "/Users/danielpalmor/CLionProjects/final_project/your_outputs/" \
                        + "spk/spk_py_input_" + str(input_num) + "_k_" + str(k[k_index]) + ".txt"
                path2 = "/Users/danielpalmor/CLionProjects/final_project/tests_inputs_outputs/outputs/spk/" \
                        + "spk_py_input_" + str(input_num) + "_k_" + str(k[k_index]) + ".txt"
            file1 = pd.read_csv(path1)
            file2 = pd.read_csv(path2)

            print("goal = " + goals[index] + ", " + c_or_py + ", file = input" + str(input_num))
            df_of_differences = file1.compare(file2, align_axis=0, keep_equal=True)
            print(df_of_differences)
            if len(df_of_differences) > 0:
                output.write("goal = " + goals[index] + ", " + c_or_py + ", file = input" + str(input_num)
                             + ", k = " + str(k[k_index]) + "\n")
                df_of_differences.replace(to_replace=np.nan, value="      ", inplace=True)
                df_of_differences.replace(to_replace=0.0, value="0.0000", inplace=True)
                df_of_differences.to_csv(output)
                output.write("\n\n")
                num_of_different_files += 1
                if input_num == 15:
                    print("\ndifferent results in: "+str(num_of_different_files)+" files")
                    output.write("\ndifferent results in: "+str(num_of_different_files)+" files")
                    output.close()
            print("")
            if input_num == 15 and k[k_index] == 5:
                print("\ndifferent results in: "+str(num_of_different_files)+" files")
                output.write("\ndifferent results in: "+str(num_of_different_files)+" files")
                output.close()

            if goals[index] == "spk":
                if twelve == input_num:
                    if mod % 3 == 0:
                        input_num += 1
                        mod = 0
                else:
                    if mod % 2 == 0:
                        input_num += 1
                k_index += 1
                mod += 1
            else:
                input_num += 1

        except FileNotFoundError:
            if not c_finished:
                input_num = 0
                c_or_py = "py"
                c_finished = True
            elif c_finished:
                input_num = 0
                c_or_py = "c"
                index += 1
                c_finished = False
        except ValueError:
            if goals[index] == "jacobi":
                output.write("goal = " + goals[index] + ", " + c_or_py + ", file = input" + str(input_num)
                             + ", k = " + str(k[k_index]) + "\n")
                output.write("our output:\n")
                output.write(file1.to_string())
                output.write("\n her output:\n")
                output.write(file2.to_string())
                output.write("\n\n")
                num_of_different_files += 1
            else:
                print("k = " + str(len(file1.columns)))
                print("files had different k\n")
                output.write("goal = " + goals[index] + ", " + c_or_py + ", file = input" + str(input_num)
                             + ", k = " + str(k[k_index]) + "\n")
                output.write("files had different sizes, our k = " + str(len(file1.columns)) + "\n\n")
                num_of_different_files += 1
            if goals[index] == "spk":
                if twelve == input_num:
                    if mod % 3 == 0:
                        input_num += 1
                        mod = 0
                else:
                    if mod % 2 == 0:
                        input_num += 1
                k_index += 1
                mod += 1
            else:
                input_num += 1
    print("\ndifferent results in: "+str(num_of_different_files)+" files")
    output.write("\ndifferent results in: "+str(num_of_different_files)+" files")
    output.close()
    output.close()


if __name__ == "__main__":
    main()
