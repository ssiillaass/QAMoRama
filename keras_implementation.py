import os
from scipy import io
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.datasets import mnist
from tensorflow.keras.utils import to_categorical
from plotting import *

def bin2dec(bin_array):
    dec_array = np.zeros([len(bin_array),1])
    for i in range(len(bin_array)):
        for j in range(len(bin_array[i])):
            power = (5-j)
            product = bin_array[i][j]*(2**power)
            dec_array[i] += product
    return dec_array

# Dataset variables
trainingsize        = 20000
testingsize         = 10000
features            = 2  #Eingang in das Netz -> [re, im]
classes             = 2 #64 verschiedene Bit_Kombinationen möglich?
hidden_layer_1_size = 1024 #
hidden_layer_2_size = 128 #
epochs              = 1000
nodes               = [features, hidden_layer_1_size,hidden_layer_2_size, classes] 

###################################################
Workspace = io.loadmat('Move_it_daten/DUL_Training_0dBm_512/2020-Apr-23_1523/dcData.mat')
#print(Workspace.keys())
MatlabCell = Workspace["normalize__2_out"]
PRMS = Workspace["prms_out"]

RxValues     = MatlabCell[0][0][0][0][0:30000] #erste 30.000 Werte komplexe Zahlen untereinander [30k, ]
TxData       = PRMS[0][0][0][0][0][0].T[0:30000] #30.000 6Bit Kombinationen untereinander (beachte Transpose) 
#Loop um die komplexen Zahlen in Real und Imag Feature zu teilen 
RxData = np.zeros([0, 2])       
for line in RxValues:                        
    RxData = np.append(RxData, [[line.real, line.imag]], axis=0)
###################################################

#Outputdaten
TxTrain            = TxData[0              :trainingsize                 ,:].astype(np.float32) 
TxTest             = TxData[trainingsize   :trainingsize+testingsize     ,:].astype(np.float32) 

TxTrain = to_categorical(bin2dec(TxTrain),num_classes=classes).astype(np.float32) 
TxTest =  to_categorical(bin2dec(TxTest) ,num_classes=classes).astype(np.float32) 

#Inputdaten
RxTrain             = RxData[0              :trainingsize                 ,:].astype(np.float32)
RxTest              = RxData[trainingsize   :trainingsize+testingsize     ,:].astype(np.float32) 

###################################################

model = tf.keras.models.Sequential([
  tf.keras.layers.Input(shape=(2,)),  
  tf.keras.layers.Dense(1024,activation='sigmoid'),
  #tf.keras.layers.Dense(256,activation='sigmoid'),
  tf.keras.layers.Dense(2, activation='softmax')
])

model.compile(
    loss='categorical_crossentropy',
    optimizer=tf.keras.optimizers.Adam(0.005),
    metrics=['accuracy'],
)

model.fit(
    RxTrain,TxTrain,
    epochs=epochs,
    validation_data=(RxTest,TxTest)
)


