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
UninstallDisplayIcon={app}\AmplifyP.exe
LicenseFile=LICENSE
SourceDir=.
OutputBaseFilename=amplifyp-windows-setup-{#Version}
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
DisableWelcomePage=no
DisableDirPage=no
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "startmenuicon"; Description: "Create a Start Menu shortcut"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "build\AmplifyP\*"; DestDir: "{app}"; Flags: recursesubdirs
Source: "src\assets\icon.ico"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\AmplifyP"; Filename: "{app}\AmplifyP.exe"; IconFilename: "{app}\icon.ico"; Tasks: startmenuicon
Name: "{autodesktop}\AmplifyP"; Filename: "{app}\AmplifyP.exe"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\AmplifyP.exe"; Description: "Launch AmplifyP"; Flags: nowait postinstall skipifsilent
