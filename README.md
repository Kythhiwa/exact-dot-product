# Точное вычисление скалярного произведения с анализом ошибок


Проект реализует алгоритм точного вычисления скалярного произведения векторов с анализом ошибок округления и сравнением с эталонными значениями (MPFR).

## 🌟 Особенности

- Точное вычисление скалярного произведения с минимизацией ошибок
- Поддержка чисел с разными порядками величин (от 1e-152 до 1e152)
- Сравнение с эталонными вычислениями с использованием MPFR
- Анализ битовых представлений результатов

## 📊 Точность достигается благодаря

1. Разбиению чисел на суммы двух с 26 знаками в мантисе.
2. Хранению промежуточных результатов с помощью 3 таблиц экспонент.

## 📦 Зависимости

- Компилятор C (gcc/clang)
- Библиотека MPFR (для эталонных вычислений)
- Стандартные математические библиотеки

## 🛠️ Сборка и запуск

1. Установите MPFR:
```bash
sudo apt-get install libmpfr-dev  # для Ubuntu/Debian
