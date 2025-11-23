import re

# Укажите имя вашего файла
filename = "output.txt"

with open(filename, 'r', encoding='utf-8') as file:
    content = file.read()

# Находим все числа с плавающей точкой после "Время выполнения: "
times = re.findall(r'Время выполнения: (\d+\.\d+) секунд', content)
times = [float(t) for t in times]

if times:
    average = sum(times) / len(times)
    print(f"Найдено измерений: {len(times)}")
    print(f"Среднее время: {average:.6f} секунд")
    print(f"Минимальное: {min(times):.6f} секунд")
    print(f"Максимальное: {max(times):.6f} секунд")
else:
    print("Данные о времени не найдены")