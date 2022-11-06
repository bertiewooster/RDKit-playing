# This only works if you call fib() in order, e.g. fib()(2), then fib()(3), fib()(4), etc.
def fib():
    x1 = 0
    x2 = 1

    def get_next_number():
        nonlocal x1, x2
        x3 = x1 + x2
        x1, x2 = x2, x3
        return x3

    return get_next_number


fibonacci = fib()

for i in range(2, 10, 2):
    num = fibonacci()
    print(f"The {i}th Fibonacci number is {num}")

# j = 7
# num = fibonacci()
# print(f"The {i}th Fibonacci number is {num}")
