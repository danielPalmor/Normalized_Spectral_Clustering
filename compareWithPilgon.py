import numpy as np
import pandas as pd


def main():
    goals = ["wam", "ddg", "lnorm", "jacobi", "spk"]
    k = [0, 3, 0, 2, 0, 3, 0, 4, 0, 10, 0, 5, 0, 20, 0, 7, 0, 10, 0, 19, 0, 25, 0, 30, 0, 10, 40, 0, 2, 0, 2, 0, 5]
    c_or_py = "c"

    k_index = 0
    input_num = 0
    index = 0
    mod = 1
    num_of_different_files = 0
    output = open("/Users/danielpalmor/Desktop/final_project/PilgonOutputs/differences.txt", "w")
    while index < len(goals):
        try:
            if goals[index] == "spk":
                c_or_py = "py"

                path1 = "/Users/danielpalmor/CLionProjects/final_project/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + "_k_" \
                        + str(k[k_index]) + ".txt"
                path2 = "/Users/danielpalmor/Desktop/final_project/PilgonOutputs/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + "_k_" \
                        + str(k[k_index]) + ".txt"
            else:
                path1 = "/Users/danielpalmor/CLionProjects/final_project/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + ".txt"
                path2 = "/Users/danielpalmor/Desktop/final_project/PilgonOutputs/your_outputs/" \
                        + goals[index] + "/" + goals[index] + "_" + c_or_py + "_input_" + str(input_num) + ".txt"

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
            if goals[index] != "spk":
                input_num += 1
            else:
                if 12 == input_num:
                    if mod % 3 == 0:
                        input_num += 1
                        mod = 0
                else:
                    if mod % 2 == 0:
                        input_num += 1
                k_index += 1
                mod += 1

            print("")

        except FileNotFoundError:
            if goals[index] != "spk":
                if c_or_py == "c":
                    c_or_py = "py"
                else:
                    c_or_py = "c"
                    index += 1
                input_num = 0
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
                if 12 == input_num:
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
        except IndexError:
            print("\ndifferent results in: " + str(num_of_different_files) + " files")
            output.write("\ndifferent results in: " + str(num_of_different_files) + " files")
            output.close()
            return


if __name__ == "__main__":
    main()
