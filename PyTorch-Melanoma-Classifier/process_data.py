import os
import cv2
import matplotlib.pyplot as plt
import numpy as np

def main():
    # 50x50 pixels
    img_size = 50

    #paths to img files
    ben_training_folder = "melanoma_cancer_dataset/train/benign/"
    mal_training_folder = "melanoma_cancer_dataset/train/malignant/"

    ben_testing_folder = "melanoma_cancer_dataset/test/benign/"
    mal_testing_folder = "melanoma_cancer_dataset/test/malignant/"

    ben_training_data = []
    mal_training_data = []
    ben_testing_data = []
    mal_testing_data = []
    
    for filename in os.listdir(ben_training_folder):
        try:
            path = ben_training_folder+filename

            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

            img = cv2.resize(img, (img_size, img_size))

            img_array = np.array(img)

            ben_training_data.append([img_array, np.array([1,0])])
        except:
            pass

    for filename in os.listdir(ben_testing_folder):
        try:
            path = ben_testing_folder+filename

            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

            img = cv2.resize(img, (img_size, img_size))

            img_array = np.array(img)

            ben_testing_data.append([img_array, np.array([1,0])])
        except:
            pass

    for filename in os.listdir(mal_training_folder):
        try:
            path = mal_training_folder+filename

            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

            img = cv2.resize(img, (img_size, img_size))

            img_array = np.array(img)

            mal_training_data.append([img_array, np.array([0,1])])
        except:
            pass

    for filename in os.listdir(mal_testing_folder):
        try:
            path = mal_testing_folder+filename

            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

            img = cv2.resize(img, (img_size, img_size))

            img_array = np.array(img)

            mal_testing_data.append([img_array, np.array([0,1])])
        except:
            pass

    ben_training_data = ben_training_data[0:len(mal_training_data)]

    print(f"Benign training count: {len(ben_training_data)}")
    print(f"Malignant training count: {len(mal_training_data)}")
    print(f"Benign testing count: {len(ben_testing_data)}")
    print(f"Malignant testing count: {len(mal_testing_data)}")

    training_data = ben_training_data + mal_training_data
    np.random.shuffle(training_data)
    np.save("melanoma_training_data.npy", np.array(training_data, dtype=object))

    testing_data = ben_testing_data + mal_testing_data
    np.random.shuffle(testing_data)
    np.save("melanoma_testing_data.npy", np.array(testing_data, dtype=object))

if __name__ == "__main__":
    main()