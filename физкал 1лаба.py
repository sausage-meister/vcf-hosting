import numpy as np
import matplotlib.pyplot as plt
import math as mth
#Здесь пишем свои замеры
time=[0, 15, 30, 45, 60, 90, 120, 150, 210, 240, 300, 360, 420, 540, 660, 780, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000, 3300, 3600, 3900, 4200, 4500, 4800]
mass=[136, 150, 155, 160, 163, 168, 180, 184, 190, 200, 207, 212, 217, 225, 231, 238, 243, 251, 258, 264, 268, 272, 275, 277, 279, 281, 282, 283, 284, 284]
#Здесь пишем свои замеры
polFunctionCF = np.polyfit(time, mass, 11)
polFunction = np.poly1d(polFunctionCF)
x = np.linspace(0,4800,90)
y = np.polyval(polFunctionCF, x)

#Нахождение массы минимального радиуса
reversedMass = mass[::-1]
for i in range(len(reversedMass)):
    if reversedMass[i]==reversedMass[i-1]:
        Minmass = reversedMass[i]
        MinmassMassive = [Minmass] * len(x)
        LastPoint = len(mass)-i

#Нахождение максимального радиуса
found = False
for i in range(len(time)):
    if i>=2 and i<=4 and (mass[i]-mass[i-1])==(mass[i-1]-mass[i-2]):
        MaxRadius = time[i]
        k=(mass[i]-mass[0])/(time[i])
        MaxRadiusMassive = []
        xmaxRad = []
        for j in range(len(time)):
            if j<=i+8:
                MaxRadiusMassive.append(time[j]*k+mass[0])
                xmaxRad.append(time[j])
        found = True
        FirstPoint = i
        break

if not found:
    MaxRadius = time[1]
    MaxRadiusMassive = []
    xmaxRad = []
    for i in range(len(time)):
        if i<= 1:
            MaxRadiusMassive.append(time[i]*((mass[1]-mass[0])/(time[1]-time[0]))+mass[0])
            xmaxRad.append(time[j])
    FirstPoint = 1



#Нахождение всех производных в точках
difFunctionCF = np.polyder(polFunctionCF)
diffunction = np.poly1d(difFunctionCF)
ydx = np.polyval(difFunctionCF, time)
#Постройка по этип производным граффиков
DifPointXTime = {}
DifPointYMass = {}
_1DifPointYMass = {}
for i in range(len(time)):
    if i>FirstPoint and i<LastPoint:
        DifPointXTime[i]=[]
        _1DifPointYMass[i]=[]
        DifPointYMass[i]=[]
        j=0
        while j<=i:
            DifPointXTime[i].append(time[j])
            _1DifPointYMass[i].append(time[j]*ydx[i]+mass[i])
            j+=1
        j-=1
        KorrectedMass=_1DifPointYMass[i][j]-_1DifPointYMass[i][0]
        for u in range(len(_1DifPointYMass[i])):
            DifPointYMass[i].append(_1DifPointYMass[i][u]-KorrectedMass)

#Постройка графиков производных по выбраным точкам(выбираем по времени оседания)
NewDifPointXTime = {}
NewDifPointYMass = {}
New_1DifPointYMass = {}
for i in range(len(time)):
    if time[i] // 60 in {75, 70, 55, 35, 20, 11, 6, 3, 0.75}: #Впиши времмя оседания выбраных точек (в минутах)
        NewDifPointXTime[i]=[]
        New_1DifPointYMass[i]=[]
        NewDifPointYMass[i]=[]
        j=0
        while j<=i:
            NewDifPointXTime[i].append(time[j])
            New_1DifPointYMass[i].append(time[j]*ydx[i]+mass[i])
            j+=1
        j-=1
        KorrectedMass=New_1DifPointYMass[i][j]-New_1DifPointYMass[i][0]
        for u in range(len(New_1DifPointYMass[i])):
            NewDifPointYMass[i].append(New_1DifPointYMass[i][u]-KorrectedMass)


#Высчитываем радиусы осевших частиц(мкм)
time = np.array(time)
radius = np.sqrt((98*10**(-3)*9*1.004*10**(-3))/(2*9.8*1000*time*(2.65-0.99821)))*10**6
ReversedRadius = radius[::-1]
ReversedRadius=ReversedRadius.tolist()
print('радиусы осевших частиц:',ReversedRadius)
#Высчитываем радиусы осевших частиц конкретно для нашей выборки
FilterReversedRadius=[]
for i in range(len(ReversedRadius)):
    #Напиши время осаждения выбранных частиц за исключением тех,
    #у которых максимальный и минимальный радиусы
    if time[i] // 60 in {55, 35, 20, 11, 6}:
        FilterReversedRadius.append(radius[i])
FilterReversedRadius = FilterReversedRadius[::-1]
FilterReversedRadius.insert(0,2.465) #Добавь сюда минимальный радиус
FilterReversedRadius.append(12.327)#Добавь сюда предпоследний с конца радиус выборки
FilterReversedRadius.append(24.654)#Добавь сюда максимальный радиус
print('радиусы осевших частиц нашей выборки:',FilterReversedRadius)

#Высчитываем радиус 
print('Масса всех радиусов вплоть до минимального:',Minmass)
print('Время Максимального радиуса:',MaxRadius)
#Высчитываем точки пересечений с одринатой(массы) для выбраных точек:
Masses=[mass[0]]
DifferenceMasses=[]
for i in range(len(time)):
    #Напиши время осаждения выбранных частиц за исключением тех,
    #у которых максимальный и минимальный радиусы
    if time[i] // 60 in {55, 35, 20, 11, 6, 3}:
        Masses.append(NewDifPointYMass[i][0])
Masses.append(Minmass)
for i in range(len(Masses)):
    if i>0:
        DifferenceMasses.append(Masses[i]-Masses[i-1])
print('Сами массы:')
for i in range(len(Masses)):
    print(Masses[i])
print('Разности масс:')
for i in range(len(DifferenceMasses)):
    print(DifferenceMasses[i])
    
#Высчитываем содержание фракции
FractionСontent = []
for i in range(len(DifferenceMasses)):
    FractionСontent.append(DifferenceMasses[i]/(Minmass-mass[0])*100)
print('Содержания фракций:')
for i in range(len(FractionСontent)):
    print(FractionСontent[i])
#Высчитываем суммарное содержание фракций
SumFractionСontent = [0.00, 19.11, 31.22, 41.04, 59.86, 70.23, 84.32, 100.00] #Впиши сюда суммарные содержания фракций, которые высчитаешь на калькуляторе

#Делим абциссу интегральной прямой на 10 равных промежутков
DifFilterReversedRadius=FilterReversedRadius[-1]-FilterReversedRadius[0]
Promejutok=DifFilterReversedRadius/10
NewFilterReversedRadius = []
for i in range(11):
    if i==0:
        NewFilterReversedRadius.append(FilterReversedRadius[i])
    else:
        NewFilterReversedRadius.append(NewFilterReversedRadius[i-1]+Promejutok)
#Высчитываем суммарное содержание фракций для новых точек на абциссе
NewSumFractionСontent = np.interp(NewFilterReversedRadius, FilterReversedRadius, SumFractionСontent)

#Высчитываем средние значения радиусов в заданном интервале (Для второй таблицы)
AverageValuesOfTheRadii = []
for i in range(len(NewFilterReversedRadius)):
    if i>0:
        AverageValuesOfTheRadii.append((NewFilterReversedRadius[i]+NewFilterReversedRadius[i-1])/2)
print('Средние значения радиуса в заданном интервалле')
for i in range(len(AverageValuesOfTheRadii)):
    print(AverageValuesOfTheRadii[i])
#Содержание фракции в данном интервале(Вторая таблица)
СontentFractionsInterval=[]
for i in range(len(NewSumFractionСontent)):
    if i>0:
        СontentFractionsInterval.append(NewSumFractionСontent[i]-NewSumFractionСontent[i-1])
print('Содержание фракций в данном интервале')
for i in range(len(СontentFractionsInterval)):
    print(СontentFractionsInterval[i])
    
#Высчитываем (m'j/r)(Вторая таблица)
MassForRadius = np.array(СontentFractionsInterval)/Promejutok
MassForRadius.tolist()
print('m`j/r')
for i in range(len(MassForRadius)):
    print(MassForRadius[i])
    
    
#Первый график
plt.figure(figsize=(8, 6))
plt.plot(time,mass,'-o', label = 'Исходные данные' ,color = 'red')
#plt.plot(x,y, label = 'Апроксимация' ,color = 'blue')
plt.plot(x,MinmassMassive, label = 'Находим массу всех фракций' ,color = 'green')
plt.plot(xmaxRad,MaxRadiusMassive, label = 'Находим максимальный радиус частиц',color = 'yellow')
plt.legend()
plt.grid(True)
plt.xlabel('T(сек)', fontsize = 14)
plt.ylabel('m(мг)', fontsize = 14)

#Второй график
plt.figure(figsize=(8, 6))
plt.plot(time,mass,'-o', label = 'Исходные данные' ,color = 'red')
for i in range(len(time)):
    if i>FirstPoint and i<LastPoint:
        plt.plot(DifPointXTime[i],DifPointYMass[i],color = 'yellow')
plt.plot(x,MinmassMassive, label = 'Находим массу каждой фракции' ,color = 'yellow')
plt.legend()
plt.grid(True)
plt.xlabel('T(сек)', fontsize = 14)
plt.ylabel('m(мг)', fontsize = 14)

#Второй график только с выбранными точками
plt.figure(figsize=(8, 6))
plt.plot(time,mass,'-o', label = 'Исходные данные' ,color = 'red')
for i in range(len(time)):
    #Напиши время осаждения выбранных частиц за исключением тех,
    #у которых максимальный и минимальный радиусы
    if time[i] // 60 in {55, 35, 20, 11, 6, 3}:
        plt.plot(NewDifPointXTime[i],NewDifPointYMass[i],color = 'yellow')
plt.plot(x,MinmassMassive, label = 'Находим массу каждой фракции' ,color = 'yellow')
plt.legend()
plt.grid(True)
plt.xlabel('T(сек)', fontsize = 14)
plt.ylabel('m(мг)', fontsize = 14)


#График интегральной кривой распределения частиц по размерам
plt.figure(figsize=(8, 6))
plt.plot(FilterReversedRadius,SumFractionСontent,'-o',color = 'red')
DifXGorizontal={}
DifYGorizontal={}
DifXVertical={}
DifYVertical={}
for i in range(len(FilterReversedRadius)):
    DifXGorizontal[i]=[]
    DifYGorizontal[i]=[]
    DifXVertical[i]=[]
    DifYVertical[i]=[]
    j=0
    while j<=i:
        DifXGorizontal[i].append(FilterReversedRadius[j])
        DifYGorizontal[i].append(SumFractionСontent[i])
        DifXVertical[i].append(FilterReversedRadius[i])
        DifYVertical[i].append(SumFractionСontent[j])
        j+=1
for i in range(len(FilterReversedRadius)):
    plt.plot(DifXGorizontal[i],DifYGorizontal[i],color = 'gray', linestyle='--')
    plt.plot(DifXVertical[i],DifYVertical[i],color = 'gray', linestyle='--')
plt.grid(True)
plt.xlabel('r(мкм)', fontsize = 14)
plt.ylabel('m`(%)', fontsize = 14)
plt.title('Интегральная кривая распределения частиц по размерам')

#График интегральной кривой распределения частиц по размерам c одинаковыми интевалами
plt.figure(figsize=(8, 6))
plt.plot(NewFilterReversedRadius,NewSumFractionСontent,'-o',color = 'red')
NewDifXGorizontal={}
NewDifYGorizontal={}
NewDifXVertical={}
NewDifYVertical={}
for i in range(len(NewFilterReversedRadius)):
    NewDifXGorizontal[i]=[]
    NewDifYGorizontal[i]=[]
    NewDifXVertical[i]=[]
    NewDifYVertical[i]=[]
    j=0
    while j<=i:
        NewDifXGorizontal[i].append(NewFilterReversedRadius[j])
        NewDifYGorizontal[i].append(NewSumFractionСontent[i])
        NewDifXVertical[i].append(NewFilterReversedRadius[i])
        NewDifYVertical[i].append(NewSumFractionСontent[j])
        j+=1
for i in range(len(NewFilterReversedRadius)):
    plt.plot(NewDifXGorizontal[i],NewDifYGorizontal[i],color = 'gray', linestyle='--')
    plt.plot(NewDifXVertical[i],NewDifYVertical[i],color = 'gray', linestyle='--')
plt.grid(True)
plt.xlabel('r(мкм)', fontsize = 14)
plt.ylabel('m`(%)', fontsize = 14)
plt.title('Интегральная кривая распределения частиц по размерам c одинаковыми интервалами')

#Строим дифференциальную кривую распределения частиц
plt.figure(figsize=(8, 6))
del NewFilterReversedRadius[-1]
plt.bar(NewFilterReversedRadius, MassForRadius,
        width=Promejutok,
        align='edge',     
        edgecolor='black',
        color='white')
plt.grid(True)
plt.xlabel('r(мкм)', fontsize = 14)
plt.ylabel('∆m`i/∆r', fontsize = 14)
plt.title('Дифференциальная кривая распределения частиц')
x_edges = np.asarray(NewFilterReversedRadius, dtype=float)
heights = np.asarray(MassForRadius, dtype=float)
x_centers = x_edges + Promejutok / 2
y_tops = heights
x_left_bottom = x_edges[0]
y_left_bottom = 0
x_right_bottom = x_edges[-1] + Promejutok
y_right_bottom = 0
x_line = np.concatenate([[x_left_bottom], x_centers, [x_right_bottom]])
y_line = np.concatenate([[y_left_bottom], y_tops,    [y_right_bottom]])
AprCF=np.polyfit(x_line, y_line, 10)
NewX_line = np.linspace(x_line[0], x_line[-1], 100)
NewY_line = np.polyval(AprCF,NewX_line)
plt.plot(NewX_line, NewY_line, color='red')

#Строим дифференциальную кривую распределения частиц в логарифмических координатах
NewXX_line = np.array(NewX_line)
NewXXX_line = np.log(NewXX_line)
NewXXX_line.tolist()

plt.figure(figsize=(8, 6))
plt.grid(True)
plt.xlabel('r(мкм)', fontsize = 14)
plt.ylabel('∆m`i/∆r', fontsize = 14)
plt.title('Дифференциальная кривая распределения частиц в логарифмических координатах')
plt.plot(NewXXX_line, NewY_line, '-o', color='red')


print('Интервалл размеров частиц отдельных фракций:',Promejutok)
plt.show()