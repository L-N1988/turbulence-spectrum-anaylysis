# Set the path to the folder containing the EPS files
$folderPath = "./concate_figures/"

# Change directory to the folder containing EPS files
Set-Location $folderPath

# Get all EPS files in the folder
$epsFiles = Get-ChildItem -Filter "pre-PSD-concate-*.eps"

# Custom sort function to sort files numerically based on the number in their names
$sortedEpsFiles = $epsFiles | Sort-Object {
    [int]([regex]::Match($_.Name, '\d+').Value)
} -Descending

# Convert each EPS file to a PNG file using ImageMagick
$pngFiles = @()
foreach ($file in $sortedEpsFiles) {
    $pngFile = [System.IO.Path]::ChangeExtension($file.FullName, ".png")
    magick convert $file.FullName $pngFile
    $pngFiles += $pngFile
}

# Create a text file containing the list of PNG files for FFmpeg
$ffmpegListFile = "file_list.txt"
$pngFiles | ForEach-Object { "'$_'" } | Out-File -FilePath $ffmpegListFile -Encoding utf8

# Use FFmpeg to create an animated GIF from the PNG files
$finalGifPath = "./prePSD-animation.gif"
ffmpeg -f concat -safe 0 -i $ffmpegListFile -vf "fps=10,scale=320:-1:flags=lanczos" -loop 0 $finalGifPath

# Clean up temporary files
# Remove-Item $ffmpegListFile
# foreach ($pngFile in $pngFiles) {
#     Remove-Item $pngFile
# }

Write-Output "Animated GIF created at $finalGifPath"
