import numpy as np
import sys
from multiprocessing import Pool
import os
import time
os.environ["CUDA_VISIBLE_DEVICES"]="0"    
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
import tensorflow as tf
import subprocess as sp
from tensorflow import keras
import subprocess as sp

def make_predictions(batch_i):
    #config = tf.compat.v1.ConfigProto()
    #config.gpu_options.per_process_gpu_memory_fraction = 0.1
    #session = , ...)
    #with tf.compat.v1.Session(config=config) as sess:
    input_fpt=np.load(fpt_root+"/pc_"+str(batch_i)+"_cicular4.npy")
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:

      # Restrict TensorFlow to only allocate 1GB of memory on the first GPU
      #try:
      tf.config.experimental.set_virtual_device_configuration(
          gpus[0],
          [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=512)])
      logical_gpus = tf.config.experimental.list_logical_devices('GPU')
        #print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
      #except RuntimeError as e:
        # Virtual devices must be set before GPUs have been initialized
      #  print(e)

    cric4_activity=keras.models.load_model("/home/jacob/More_Data/ligand_NN_2_15_21/cicular_4_models_6_17_21/model_"+str(fpt_type)+".tf",custom_objects={'crossentropy':crossentropy})
    cric4_activity.compile(
            loss=keras.losses.MeanSquaredError(),#keras.losses.MeanSquaredError(),
            optimizer=keras.optimizers.Adam(learning_rate=1),
            metrics=["MeanSquaredError"],
        )
    bsize=min([input_fpt.shape[0],1000])
    #print(input_fpt.shape)
    if input_fpt.shape[0] !=0:
        pred_act=cric4_activity.predict(input_fpt,batch_size=bsize)
        #load smiles 
        input_smi=np.load(fpt_root+"/pc_"+str(batch_i)+"_smiles.npy")
        #np.save("/home/jacob/More_Data/ligand_NN_2_15_21/FDB-17/for_paper/fdb_17_"+str(batch_i)+"_pred",pred_act)
        active_smi=input_smi[pred_act[:,0]>0.9]
        np.save(fpt_root+"/active_smi_"+str(batch_i),active_smi)
        os.remove(fpt_root+"/pc_"+str(batch_i)+"_smiles.npy")
        os.remove(fpt_root+"/pc_"+str(batch_i)+"_cicular4.npy")
        return active_smi.shape[0]
    else:
        os.remove(fpt_root+"/pc_"+str(batch_i)+"_smiles.npy")
        os.remove(fpt_root+"/pc_"+str(batch_i)+"_cicular4.npy")
        return 0
    
def crossentropy(y_true,y_pred):
    return tf.reduce_mean(-1.*y_true[:,0]*tf.math.log(tf.clip_by_value(y_pred[:,0],1e-10,1.)) - (y_true[:,1])*tf.math.log(tf.clip_by_value(y_pred[:,1],1e-10,1.)))


delta=500000

ind_i=1
n_proc=24
ind_f=delta
#ind_i=131500001
#ind_f=132000000
ind_i=1
ind_f=500000

for i in range(312):#range(312):    
    ti=time.time()
    #print(i)
    sp.call("curl -p "+"""--insecure "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_"""+f'{ind_i:09}'+"""_"""+f'{ind_f:09}'+""".sdf.gz" """+"-o Compound_"+f'{ind_i:09}'+"""_"""+f'{ind_f:09}'+".sdf.gz",shell=True)
    sp.call("gunzip Compound_"+f'{ind_i:09}'+"""_"""+f'{ind_f:09}'+".sdf.gz",shell=True)
    #print(time()-ti)
    #exit()

    fpt_root="/home/jacob/More_Data/ligand_NN_2_15_21/screening_compounds/ml_model/pubchem/pub_chem_"+str(i)
    sdf_file="Compound_"+f'{ind_i:09}'+"""_"""+f'{ind_f:09}'+".sdf"
    if not os.path.isdir(fpt_root):
        os.mkdir(fpt_root)
    #multiprocess parsing the .sdf files to circular features length 4 
    #compute fingerprints for molecules in each .sdf
    #multiprocess generating fingerprints
    if os.path.isfile("completed.txt"):
        os.remove("completed.txt") 
    fo=open("make_fpts.sh",'w')
    fo.write("""eval "$(conda shell.bash hook)" """+"\n")
    fo.write("conda activate my-rdkit-env"+"\n")
    fo.write("python generate_circ4.py "+sdf_file+" "+fpt_root+" \n")
    fo.close()
    sp.call("bash make_fpts.sh",shell=True)

    waiting=1
    while waiting:
        if os.path.isfile("completed.txt"):
            waiting=0
        else:
            time.sleep(.01)
    os.remove("completed.txt") 

    #determine n_batch
    batch_files=os.listdir(fpt_root+"/")
    batched_files=[int(x.split("/")[-1].split("_")[1]) for x in batch_files if "cicular4" in x]
    if len(batched_files)>0:
        n_batch=max(batched_files)+1
        #fpt_file_root+"/pc_"+str(i)+"_cicular4",fpt_file_root+"/pc_"+str(i)+"_smiles"
    
        #predict activity
        
        #load the tensorflow model
        fpt_type="cicular4"
    
        ##parallel predictions
        if __name__ == '__main__':
            with Pool(12) as p:
                n_active=p.map(make_predictions, range(n_batch))
    
        #save smiles if active

    os.remove(sdf_file)

    ind_i+=delta
    ind_f+=delta
    if len(batched_files)>0:
        print(i,time.time()-ti,sum(n_active))
