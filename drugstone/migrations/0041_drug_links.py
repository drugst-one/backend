# Generated by Django 3.0.5 on 2020-05-31 13:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0040_auto_20200528_1539'),
    ]

    operations = [
        migrations.AddField(
            model_name='drug',
            name='links',
            field=models.CharField(default='', max_length=1024),
        ),
    ]
