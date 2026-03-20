import random
dataset_list = ['PXD000554', 'PXD012395', 'PXD013271', 'PXD019647', 'PXD021109', 'PXD028028', 'PXD035029', 'PXD037381']

random.seed(1234)
seven_sample = random.sample(dataset_list, 3)
print(seven_sample)

dataset_list_7_1 = [i for i in dataset_list if i != seven_sample[0]]
print(dataset_list_7_1)
dataset_list_7_2 = [i for i in dataset_list if i != seven_sample[1]]
print(dataset_list_7_2)
dataset_list_7_3 = [i for i in dataset_list if i != seven_sample[2]]
print(dataset_list_7_3)
print(dataset_list)

six_sample = random.sample(dataset_list, 6)
print(six_sample)
dataset_list_6_1 = [i for i in dataset_list if i not in six_sample[0:2]]
print(dataset_list_6_1)
dataset_list_6_2 = [i for i in dataset_list if i not in six_sample[2:4]]
print(dataset_list_6_2)
dataset_list_6_3 = [i for i in dataset_list if i not in six_sample[4:6]]
print(dataset_list_6_3)

random.shuffle(dataset_list)
print(dataset_list)
dataset_list_5_1 = dataset_list[0:5]
print(dataset_list_5_1)
dataset_list_5_2 = dataset_list[1:3] + dataset_list[5:8]
print(dataset_list_5_2)
dataset_list_5_3 = dataset_list[3:8]
print(dataset_list_5_3)
