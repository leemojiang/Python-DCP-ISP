
import rawpy
import numpy as np
import matplotlib.pyplot as plt
import time

def gamma(x,colorspace='sRGB'): #Gamma变换
    '''
        x
    '''


    y=np. zeros (x. shape)
    y[x>1]=1
    if colorspace in ( 'sRGB', 'srgb'):
        y[(x>=0)&(x<=0.0031308)]=(323/25*x[ (x>=0)&(x<=0.0031308)])
        y[(x<=1)&(x>0.0031308)]=(1.055*abs(x[ (x<=1)&(x>0.0031308)])**(1/2.4)-0.055)
    elif colorspace in ('my'):  
        y[ (x>=0)&(x<=1)]=(1.42*(1-(0.42/(x[(x>=0)&(x<=1)]+0.42))))
    elif colorspace in ('P3'):  
        y[ (x>=0)&(x<=1)]=x[ (x>=0)&(x<=1)]**(1/2.6)
    elif (type(colorspace)==float)|(type(colorspace)==int):
        beta=colorspace
        y[ (x>=0)&(x<=1)]=((1+beta)*(1-(beta/(x[(x>=0)&(x<=1)]+beta))))
    return y

def gamma_reverse(x,colorspace= 'sRGB'): #逆Gamma变换
    y=np.zeros(x.shape)
    y[x>1]=1
    if colorspace in ('sRGB', 'srgb'):
        y[(x>=0)&(x<=0.04045)]=x[(x>=0)&(x<=0.04045)]/12.92
        y[(x>0.04045)&(x<=1)]=((x[(x>0.04045)&(x<=1)]+0.055)/1.055)**2.4
    elif colorspace in ('my'):
        y[(x>=0)&(x<=1)]=0.42/(1-(x[(x>=0)&(x<=1)]/1.42))-0.42         
    return y


def awb(img, awb_para, normalize = True):  #图像做白平衡
    '''
        img [H,W,3] dim2 in R G B channal
        awb_para [R_gain,G_gain,B_gain]

        output: white balanced image in shape [H,W,3] [0,1+]
        
        Don't use normalize
    '''
    img_ = img.copy()

    img_[...,0] =  img[...,0]*awb_para[0] # R_gain
    img_[...,1] =  img[...,1]*awb_para[1] # G_gain
    img_[...,2] =  img[...,2]*awb_para[2] # B_gain

    if normalize:
        print("Only use normalize For testing!")
        img_ = img_/img_.max()

    return img_

def read_raw_simple(file_path):
    # 如果用 with rawpy.imread('your_file.arw') as H_raw: 的写法 出了with block 会报错
    H_raw=rawpy.imread(file_path)

    # 读取raw data
    raw=H_raw.raw_image
    print(f"RAW Shape:{raw.shape}")

    OB=H_raw.black_level_per_channel[0]
    # white_level=H_raw.camera_white_level_per_channel[0]
    white_level = H_raw.white_level
    print(f"OB {OB} \n white_level {white_level} \n MAX{raw.max()}")

    img=bayer_demosaic(raw)
    img[img<OB]=OB
    img=(img-OB).astype('float32')/(white_level-OB)

    img[img < 0] = 0

    return img,H_raw

def bayer_demosaic(raw,bayer='RG'): #朴素的bayer插值算法
    if bayer=='RG':
        img_r=raw[0::2,0::2]
        img_gr=raw[0::2,1::2]
        img_gb=raw[1::2,0::2]
        img_b=raw[1::2,1::2]
    img=np.dstack((img_r,(img_gr+img_gb)/2,img_b))
    return img   

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

def im2vector(img): #将图片转换为向量形式
    '''
        img [H,W,3]
        rgb [H*W,3]
        func_reverse reshape rgb to img
    '''
    size=img.shape
    rgb=np.reshape(img,(size[0]*size[1],3))
    func_reverse=lambda rgb : np.reshape(rgb,(size[0],size[1],size[2]))
    return rgb, func_reverse    

def manual_wb_area(img):
    '''
        img [H,W,3]

        list of selected wb area ([xmin,xmax,ymin,ymax],rgb_mean)
    '''
    plt.figure()
    plt.imshow(img)
    
    tellme("Select 2 points for a region to calculate WB, click to begin")
    plt.waitforbuttonpress()

    wb_list = []

    while True:
        pts = []
        while len(pts) < 2:
            tellme('Select 2 corners with mouse')
            pts = np.asarray(plt.ginput(2, timeout=-1))
            if len(pts) < 2:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second
        
        xmin = int(min(pts[:,0]))
        ymin = int(min(pts[:,1]))
        xmax = int(max(pts[:,0]))
        ymax = int(max(pts[:,1]))
        
        ph = plt.fill([xmin,xmax,xmax,xmin], [ymin,ymin,ymax,ymax], 'r', lw=2,alpha=0.3)
        
        block = img[ymin:ymax,xmin:xmax,:] #[H1,W1,3]
        rgb_vector,_=im2vector(block) #[H1*W1,3]
        rgb_mean=rgb_vector.mean(axis=0) #[3]
        rgb_std=rgb_vector.std(axis=0) #[3]
        
        print(f"Block at X:{xmin} {xmax} Y:{ymin} {ymax} Total Pixel {rgb_vector.shape[0]}")
        print(f"Mean {rgb_mean}")
        print(f"STD {rgb_std}")

        wb_list.append(([xmin,xmax,ymin,ymax],rgb_mean))

        tellme('Key click for quit, mouse click for extra data')
        if plt.waitforbuttonpress():
            tellme('Finish')
            break

        # # Get rid of fill
        # for p in ph:
        #     p.remove()
    tellme('Done')
    return wb_list