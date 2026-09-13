[Setup]
AppName=AmplifyP
AppVersion={#Version}
AppPublisher=Fufu Fang
AppPublisherURL=https://github.com/fangfufu/AmplifyP
AppSupportURL=https://github.com/fangfufu/AmplifyP/issues
AppUpdatesURL=https://github.com/fangfufu/AmplifyP/releases
DefaultDirName={localappdata}\Programs\AmplifyP
DefaultGroupName=AmplifyP
SetupIconFile=src/assets/icon.ico
UninstallDisplayIcon={app}\icon.ico
LicenseFile=LICENSE
SourceDir=.
OutputBaseFilename=amplifyp-windows-setup-{#Version}
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
DisableWelcomePage=no
DisableDirPage=no
AllowNoIcons=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[InstallDelete]
Type: filesandordirs; Name: "{app}\DLLs"
Type: filesandordirs; Name: "{app}\Lib"
Type: filesandordirs; Name: "{app}\app"
Type: filesandordirs; Name: "{app}\data"
Type: filesandordirs; Name: "{app}\site-packages"
Type: files; Name: "{app}\*.dll"
Type: files; Name: "{app}\*.pyd"

[Files]
Source: "build\AmplifyP\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "src\assets\icon.ico"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\AmplifyP"; Filename: "{app}\AmplifyP.exe"; IconFilename: "{app}\icon.ico"; WorkingDir: "{app}"
Name: "{autodesktop}\AmplifyP"; Filename: "{app}\AmplifyP.exe"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon; WorkingDir: "{app}"

[Run]
Filename: "{app}\AmplifyP.exe"; Description: "Launch AmplifyP"; Flags: nowait postinstall skipifsilent; WorkingDir: "{app}"
