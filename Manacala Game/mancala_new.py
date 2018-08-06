# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 21:27:28 2018

@author: ahmed mubin Math Finance
"""

import numpy as np

mylist=[0,4,4,4,4,4,4,0,4,4,4,4,4,4]

idx=[0,1,2,3,4,5,6,7,13,12,11,10,9,8]

macalaB_idx=0
macalaA_idx=7

def create_board():
    print("****","Board Layout Start","****")
    print("----------------------------------------------------------")
    print(" cup number player_2->","12,11,10, 9, 8, 7")
    print("************************************************")
    print("                        ",str(mylist[8:14]).strip('[]'))
    print("                      ",mylist[0],"                  ",mylist[7])
    print("                        ",str(mylist[1:7]).strip('[]'))
    print("**************************************************")
    print(" cup number player_1->"," 1, 2, 3, 4, 5, 6")
    print("----------------------------------------------------------")
    print("****","Board Layout end","****")

def start_play():
    create_board()
    player1_play()

def player1_play():
    cup_num=input("Player 1 Enter the cup number=")
    while(checkPlayer1Input(int(cup_num))!=True):
        cup_num=input("Player 1 Enter the correct cup number marked between 1 and 6 on console=")
    player1_move(int(cup_num))

def player2_play():
    cup_num=input("Player 2 Enter the cup number=")
    while(checkPlayer2Input(int(cup_num))!=True):
        cup_num=input("Player 2 Enter the correct cup number marked between 7 and 12 on console=")
    player2_move(int(cup_num))

def checkPlayer1Input(cupnum):
    if(cupnum >=1) and (cupnum <=6):
        return True
    else:
        return False

def checkPlayer2Input(cupnum):
    if(cupnum >=7) and (cupnum <=12):
        return True
    else:
        return False

def player1_move_again():
    #checkwin_plA()
    player1_play()

def reset_idx(j,cnt):
    flagA=0
    temp=cnt;
    print("reset_idx",j," ",cnt)
    for i in range(j,j+cnt):
        temp=temp-1
        if(i<=13):
            mylist[idx[i]]=mylist[idx[i]]+1
        else:
            if(mylist[idx[i-13]]==0 & temp==0):
                flagA=1                
            mylist[idx[i-13]]=mylist[idx[i-13]]+1
            if (mylist[idx[i-6]] > 0) and (flagA==1):
                mylist[macalaA_idx]=mylist[macalaA_idx]+mylist[idx[i-6]]+mylist[idx[i-13]]
                mylist[idx[i-6]]=0
                mylist[idx[i-13]]=0
                flagA=0
        
   # if(i == macalaA_idx):
  #      player1_move_again()  

def player1_move(cupIdx):
    tempFlag=0
    cnt=mylist[cupIdx]
    if(cnt!=0):
        mylist[cupIdx]=0
        for i in range(cupIdx+1,cupIdx+cnt+1):          
            if(i<=13):
                if(checkifcupEmpty(mylist[i])==True) and (i!=macalaA_idx):
                    tempFlag=1
                mylist[i]=mylist[i]+1
                if tempFlag == 1:
                    tempFlag=0
                    if(mylist[i+7]>0) and(cnt <=1):
                        mylist[macalaA_idx]=mylist[macalaA_idx]+mylist[i]+mylist[i+7]
                        mylist[i]=0
                        mylist[i+7]=0             
            else:
                mylist[i-13]=mylist[i-13]+1
            cnt=cnt-1
            if (i == macalaA_idx) and (cnt>=1):
                reset_idx(i+1,cnt)
                break;
        create_board()
        if(isGameOver() == True):
            _printResult()
            return
        else:
            if(i == macalaA_idx) and (cnt==0):
                player1_move_again()
            else:
                player2_play()
    else:
        print("Select the cup with marbles in it, not the empty cups")
        player1_play()




def _printResult():
    if(mylist[macalaB_idx] > mylist[macalaA_idx]):
        print("Player 2 won by ",mylist[macalaB_idx],"to ",mylist[macalaA_idx])
    else:
        print("Player 1 won by ",mylist[macalaA_idx],"to ",mylist[macalaB_idx])
    


def isGameOver():
    if(Is_Player1GameOver()==True):
        player2_Addthemarbes()
    elif(Is_Player2GameOver()==True): 
        player1_Addthemarbes()
    else:
        return False    
    return True
        

        
def player2_Addthemarbes():
    for i in range(8,14):
        mylist[macalaB_idx] = mylist[macalaB_idx] + mylist[i]
        

def player1_Addthemarbes():
    for i in range(1,7):
        mylist[macalaA_idx] = mylist[macalaA_idx] + mylist[i]   
        
def Is_Player1GameOver():
    for i in range(1,7):
        if mylist[i] == 0:
            continue
        else:
            return False
    return True

def Is_Player2GameOver():
    for i in range(8,14):
        if mylist[i] == 0:
            continue
        else:
            return False
    return True
    
def reset_idx_pl21(j,cnt):
    flagB=0
    temp=cnt;
    print("reset_idx",j," ",cnt)
    for i in range(j,j+cnt):
        if(i<=13):
            if(idx[i]!=macalaA_idx):
                if (flagB==1) and (mylist[idx[i]]==0):
                    if(mylist[idx[i]-7]>0):
                        mylist[macalaB_idx]=mylist[macalaB_idx]+mylist[idx[i]]+mylist[idx[i-7]]
                        mylist[idx[i-6]]=0
                        mylist[idx[i-13]]=0
                        flagB=0
                else:
                    mylist[idx[i]]=mylist[idx[i]]+1
            else:
                flagB=1
        temp=temp-1
        
def reset_idx_pl2(j,cnt):
    flag=0
    temp=cnt
    print("reset_idx",j," ",cnt)
    for i in range(j,j+cnt):
        temp=temp-1
        if(i<=13):
            if (idx[i] > 6):
                if(checkifcupEmpty(mylist[idx[i+1]])==True) and (temp==0):
                    flag=1
                mylist[idx[i+1]]=mylist[idx[i+1]]+1
                if(flag == 1):
                    if(mylist[idx[i+1]-7]>0):
                        mylist[macalaB_idx]=mylist[macalaB_idx]+mylist[idx[i+1]]+mylist[idx[i+1]-7]
                        mylist[idx[i+1]]=0
                        mylist[idx[i+1]-7]=0
            else:
                mylist[idx[i]]=mylist[idx[i]]+1
                
  #  if(i == macalaA_idx):
    #    player1_move_again() 

def player2_move_again():
    player2_play()

def checkifcupEmpty(num):
    if num == 0:
        return True
    else:
        return False
    
        
def player2_move(cupIdx):
    tempFlag=0
    cnt=mylist[idx[cupIdx+1]]
    if(cnt!=0):
        mylist[idx[cupIdx+1]]=0
        for i in range(cupIdx+2,cupIdx+cnt+2):
            if(i<=13):
                if(checkifcupEmpty(mylist[idx[i]]) == True) and (i!=macalaB_idx):
                    tempFlag = 1
                mylist[idx[i]]=mylist[idx[i]]+1
                if tempFlag==1:
                    tempFlag=0
                    if(mylist[idx[i]-7]>0) and (cnt <=1) :
                        mylist[macalaB_idx]=mylist[macalaB_idx]+mylist[idx[i]]+mylist[idx[i]-7]
                        mylist[idx[i]]=0
                        mylist[idx[i]-7]=0
            else:
                mylist[idx[i-14]]=mylist[idx[i-14]]+1
            cnt=cnt-1
            if (idx[i-14]==macalaB_idx) and (cnt>=1):
                reset_idx_pl2(idx[i-14]+1,cnt)
                break            
        create_board()
        if(isGameOver() == True):
            _printResult()
            return        
#        print("debug",i)
        if(idx[i-14]==macalaB_idx) and (cnt==0):
            player2_move_again()
        else:
            #checkwin_plA()
#            print("calling player1_move 2" )
#            cup_num=input("Player 1 Enter the cup number=")
#            player1_move(int(cup_num))   
            player1_play()    
    else:
        print("Select the cup with marbles in it, not the empty cups")
        player2_play()

start_play()