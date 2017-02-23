def merge_sort(arr):
    n = len(arr)
    # Base case
    if n <= 1:
        return arr

    # Recursive case
    split_idx = int(n / 2) # Round down for odd n

    # Recursively sort the left and right halves of the array
    left_sorted = merge_sort(arr[:split_idx])
    right_sorted = merge_sort(arr[split_idx:])

    # Join the sorted halves in order
    res = []
    i, j = 0, 0
    while i < len(left_sorted) and j < len(right_sorted):
        if left_sorted[i] <= right_sorted[j]:
            res.append(left_sorted[i])
            i += 1
        else:
            res.append(right_sorted[j])
            j += 1
    while i < len(left_sorted):
        res.append(left_sorted[i])
        i += 1
    while j < len(right_sorted):
        res.append(right_sorted[j])
        j += 1
    return res
